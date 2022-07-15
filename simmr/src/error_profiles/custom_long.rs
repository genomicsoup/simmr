/**
 * file: custom_long.rs
 * desc: A custom long-read error profile. This generates phred scores and errors using
 *       an empirical distribution constructed from observations present in a distribution file.
 *       Read lengths are generated using a gamma distribution.
 *       
 */
use rand::prelude::SliceRandom;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Distribution, Gamma, Normal, WeightedAliasIndex};
use stats::stddev;
use std::collections::HashMap;
use std::f64::consts::PI;

use shared::encoding;

use super::base;
use crate::util;

pub struct CustomLongErrorProfile {
    // These are NOT unused by this profile
    //pub read_length: u16,
    //pub insert_size: u16,
    pub bin_size: u8,
    pub quality_bins: Vec<Vec<u32>>,
    pub quality_bins2: Vec<Vec<f32>>,
    pub kmer_size: usize,
    pub kmer_probabilities: HashMap<u32, Vec<(u32, f32)>>,
}

// TODO: put this elsewhere
// Use the normal kernel to generate a value from the density function using the given value (x).
// Implements PDF formula here: https://en.wikipedia.org/wiki/Gaussian_function
fn gaussian(x: f64, xs: &[f64], bandwidth: f64) -> f64 {
    let sum: f64 = xs
        .iter()
        .map(|f| (x - f) / bandwidth)
        .map(|f| (-0.5 * f.powi(2)).exp())
        .sum();

    sum / ((2.0 * PI).sqrt() * xs.len() as f64 * bandwidth)
}

impl base::ErrorProfile for CustomLongErrorProfile {
    /**
     * Unusable functions. These can't and shouldn't be used for long-read error profiles.
     * Calling these when simulating long reads (which should never happen) will cause a panic.
     */
    fn get_read_length(&self) -> u16 {
        panic!("get_read_length() is not usable when simulating long reads");
    }

    fn get_insert_size(&self) -> u16 {
        panic!("get_insert_size() is not usable when simulating long reads");
    }

    /**
     * Generates a random read length using a gamma distribution with a mean of 20K and
     * standard deviation of 15K.
     */
    fn get_random_read_length(&self, seed: Option<u64>) -> u16 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        // TODO: make these configurable?
        let mean: f32 = 20_000.0;
        let std_dev: f32 = 15_000.0;

        // Convert mean and standard deviation to shape and scale parameters. It's easier to
        // reason and visualize distributions when talking about mean/stddev. I have no fucking
        // clue what the shape/scale of a distribution means...
        let shape = (mean / std_dev).powf(2.0);
        let scale = std_dev.powf(2.0) / mean;

        // Sample a new value, truncate the resulting float into a u16
        rng.sample(&Gamma::new(shape, scale).unwrap()).floor() as u16
    }

    /**
     * Simulates phred scores using density estimation. Each position in the simulated read
     * uses it's own density function to simulate scores at that position. If the simulated
     * read is longer than what we have positions for, we just use the last position available.
     */
    fn simulate_phred_scores(&self, seq_length: usize, seed: Option<u64>) -> Vec<u8> {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let mut scores = Vec::new();

        // Simulate quality score at each position in the read
        for i in 0..seq_length {
            // Get bins for this position
            let bins = if i > self.quality_bins.len() {
                &self.quality_bins[self.quality_bins.len() - 1]
            } else {
                &self.quality_bins[i]
            };

            let mut centers = Vec::new();

            // Since we're not using individual data points, calculate centers for each bin
            for (bndx, count) in bins.iter().enumerate() {
                // We store center * count values, since that's the amount of observations in that bin
                for _ in 0..*count {
                    // We generate random values within this bin. We could also use bin centers
                    // but I think this generates more realistic scores
                    let val = if bndx == 0 {
                        rng.gen_range(0..self.bin_size) as u32
                    } else {
                        rng.gen_range(
                            ((bndx - 1) * self.bin_size as usize)..(bndx * self.bin_size as usize),
                        ) as u32
                    };
                    centers.push(val as f64);

                    //centers.push(if bndx == 0 {
                    //    self.bin_size as f64 / 2.0
                    //} else {
                    //    (((bndx - 1) * self.bin_size as usize) + (bndx * self.bin_size as usize))
                    //        as f64
                    //        / 2.0
                    //});
                }
            }

            // Estimate the bandwidth parameter using Silverman's rule of thumb, not the best but it's fast
            let bandwidth = 1.06
                * stddev(centers.clone().into_iter())
                * (centers.len() as f64).powf(-1.0 / 5.0);

            // Generate a random quality score that corresponds to the distribution at this position
            // TODO: don't hardcode these, a better thing would be to select the min/max quality scores
            // for this position
            // TODO: also, should probably be simulating error/accuracy and then converting that to
            // a phred score...
            let qs = centers.choose(&mut rng).unwrap();
            let normal = Normal::new(*qs, bandwidth).unwrap();
            let sim_qs = normal.sample(&mut rng);

            scores.push(sim_qs.round() as u8);
            continue;

            // This is the harder way to sample, I think this can just be replaced by the above
            // code
            loop {
                let qs = rng.gen_range(0..70) as u8;

                if rng.gen::<f64>() <= gaussian(qs as f64, &centers, bandwidth) {
                    scores.push(qs as u8);
                    break;
                }
            }
        }

        scores
    }

    // Introduce errors into the given sequence by using sequencing errors observed across
    // kmers from sequencing experiments. This function kmerizes the given sequence, then for
    // every kmer, uses a weighted random distribution (which essentially serves as an
    // empirical distribution) to select an alternate kmer to replace it with. The alternate
    // kmer could also be the original kmer itself.
    fn simulate_errors(&self, sequence: &[u8], seed: Option<u64>) -> Vec<u8> {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Make a copy of the sequence which we'll mutate and introduce errors into
        let mut new_sequence = sequence.to_vec();

        for i in 0..sequence.len() {
            if (i + self.kmer_size) > sequence.len() {
                break;
            }

            // Grab the kmer
            let kmer = &new_sequence[i..i + self.kmer_size];
            // Encode it. If it can't be encoded b/c of characters other than ACGT, skip it
            let encoded_kmer = match encoding::three_bit_encode_kmer(kmer, self.kmer_size) {
                Ok(k) => k,
                Err(_) => continue,
            };

            // Find the list of alternate kmers, if this kmer hasn't been observed then move on
            let alts = match self.kmer_probabilities.get(&encoded_kmer) {
                Some(vs) => vs,
                None => continue,
            };

            // Weighted distribution sampling
            let weighted_dist =
                WeightedAliasIndex::new(alts.iter().map(|(_, w)| *w).collect()).unwrap();

            // Get an alternate kmer which could include the og kmer itself
            let alt_encoded_kmer = alts[weighted_dist.sample(&mut rng)].0;

            // Decode the kmer into a byte array
            let alt_kmer =
                encoding::three_bit_decode_kmer(alt_encoded_kmer, self.kmer_size, true).unwrap();

            // Replace this kmer with the alt in the sequence. Splice is neat in that if alt_kmer has deletions,
            // it will delete the elements at those indices from new_sequence too.
            new_sequence.splice(i..i + self.kmer_size, alt_kmer);
        }

        new_sequence
    }

    /**
     * This doesn't actually do anything for custom profiles but since it's used by other long
     * read profiles, it still needs to return the sequence.
     */
    fn simulate_point_mutations(
        &self,
        sequence: &[u8],
        _quality: &Vec<u8>,
        _seed: Option<u64>,
    ) -> Vec<u8> {
        sequence.to_vec()
    }
}

#[cfg(test)]
mod tests {

    use base::ErrorProfile;

    use super::*;

    #[test]
    fn test_gaussian_function() {
        assert!(
            gaussian(
                4.0,
                &vec![9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
                0.1
            ) == 0.3989422804014327
        );
    }

    /* Not really a test, just used to eyeball.
    #[test]
    fn test_simulate_phred_scores() {
        let bins = vec![
            vec![10, 20, 30, 40, 50, 60, 60, 40, 20, 10],
            vec![10, 20, 30, 40, 50, 60, 60, 50, 30, 20],
            vec![5, 10, 30, 40, 60, 60, 60, 60, 40, 20],
            vec![20, 20, 30, 40, 40, 30, 20, 20, 5, 5],
            vec![30, 30, 30, 30, 20, 30, 20, 20, 5, 5],
        ];
        let binsf: Vec<Vec<f32>> = bins
            .iter()
            .map(|vs| vs.iter().map(|v| *v as f32).collect())
            .collect();
        let cel = CustomLongErrorProfile {
            bin_size: 5,
            quality_bins: bins,
            quality_bins2: binsf,
            kmer_size: 7,
            kmer_probabilities: HashMap::new(),
        };
        println!("{:?}", cel.simulate_phred_scores(5, None));
    }
    */

    #[test]
    fn test_simulate_errors() {
        let cel = CustomLongErrorProfile {
            bin_size: 0,
            quality_bins: vec![],
            quality_bins2: vec![],
            kmer_size: 3,
            kmer_probabilities: HashMap::from([
                (
                    encoding::three_bit_encode_kmer("ACC".as_bytes(), 3).unwrap(),
                    vec![
                        (
                            encoding::three_bit_encode_kmer("ACC".as_bytes(), 3).unwrap(),
                            0.0,
                        ),
                        (
                            encoding::three_bit_encode_kmer("CAT".as_bytes(), 3).unwrap(),
                            1.0,
                        ),
                    ],
                ),
                (
                    encoding::three_bit_encode_kmer("ATC".as_bytes(), 3).unwrap(),
                    vec![(
                        encoding::three_bit_encode_kmer("ATC".as_bytes(), 3).unwrap(),
                        1.0,
                    )],
                ),
                (
                    encoding::three_bit_encode_kmer("TCG".as_bytes(), 3).unwrap(),
                    vec![(
                        encoding::three_bit_encode_kmer("TGT".as_bytes(), 3).unwrap(),
                        1.0,
                    )],
                ),
            ]),
        };

        //ACCCG
        //CATCG
        //CATGT
        assert!(
            std::str::from_utf8(&cel.simulate_errors("ACCCG".as_bytes(), None)).unwrap() == "CATGT"
        );
    }
}
