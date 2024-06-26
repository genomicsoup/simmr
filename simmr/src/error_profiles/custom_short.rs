/**
 * file: custom_short.rs
 * desc: A custom short-read error profile. This generates phred scores and errors using
 *       an empirical distribution constructed from observations present in a distribution file.
 *       Read lengths are generated using a gamma distribution.
 *       
 */
use rand::prelude::SliceRandom;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Distribution, Gamma, Normal, Uniform, WeightedAliasIndex};
use stats::stddev;
use std::collections::HashMap;
use std::f64::consts::PI;

use shared::encoding;

use super::base;
use crate::util;

pub struct QualityScoreDistribution {
    pub bandwidth: f64,
    pub centers: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct CustomPDF {
    // PDF density
    pub density: Vec<WeightedAliasIndex<f64>>,
    // The bin ranges for the PDF
    //bin_ranges: Vec<(u32, u32)>,
    pub per_bin_values: Vec<Vec<Uniform<u32>>>,
}

#[derive(Clone, Debug)]
pub struct CustomShortErrorProfile {
    // These are NOT unused by this profile
    //pub read_length: u16,
    //pub insert_size: u16,
    //pub bin_size: u8,
    //pub quality_bins: Vec<Vec<u32>>,
    //pub quality_bins2: Vec<Vec<f32>>,
    //pub kmer_size: usize,
    //pub kmer_probabilities: HashMap<u32, Vec<(u32, f32)>>,
    //pub insert_size_mean: f64,
    //pub insert_size_std: f64,
    //pub read_length_mean: f64,
    //pub read_length_std: f64,
    pub model_params: encoding::ErrorModelParams,
    //pub quality_score_distributions: Vec<QualityScoreDistribution>,
    pub quality_score_pdf: CustomPDF,
    pub read_length_pdf: CustomPDF,
    pub insert_size_pdf: Option<CustomPDF>,
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

impl CustomPDF {
    /**
     * Construct a custom PDF from a model file. Generates the weight rand distributions
     * needed to sample from the PDF.
     */
    pub fn new(binned_density: &Vec<encoding::Bins>) -> CustomPDF {
        // Build the PDFs for for each base pair quality score
        let density: Vec<WeightedAliasIndex<f64>> = binned_density
            .iter()
            .map(|bins| {
                let dist = WeightedAliasIndex::new(bins.binned_density.clone()).unwrap();
                dist
            })
            .collect();

        // Build a set of uniform distributions per bin
        let per_bin_values: Vec<Vec<Uniform<u32>>> = binned_density
            .iter()
            .map(|bins| {
                bins.bin_ranges
                    .iter()
                    .map(|(start, end)| Uniform::new_inclusive(*start, *end))
                    .collect()
            })
            .collect();

        Self {
            density,
            per_bin_values,
        }
    }

    /**
     * Sample a new data point from the PDF using the given index.
     * This is mainly designed for use with a distribution of quality scores, where each
     * per base quality score has its own PDF and the index corresponds to the bp position.
     * This will panic if the index is out of bounds.
     */
    pub fn sample_with_index(&self, index: usize, seed: Option<u64>) -> u32 {
        if index >= self.density.len() || index >= self.per_bin_values.len() {
            panic!("PDF index {} out of bounds", index);
        }

        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Get a pdf to sample from
        let pdf = &self.density[index];

        // Select a bin to sample from
        let sampled_bins = &self.per_bin_values[index];
        let sampled_bin = sampled_bins[pdf.sample(&mut rng)];

        // Generate a random data point using this PDF
        let a = sampled_bin.sample(&mut rng);
        //println!("Sampled {} from bin {}", a, index);
        a
    }

    /**
     * Sample a new data point from the PDF. Assumes there is only a single PDF in the density
     * vector. Mainly designed for use with anything that doesn't require multiple PDFs, e.g.,
     * read lengths and insert sizes.
     */
    pub fn sample(&self, seed: Option<u64>) -> u32 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Get a pdf to sample from
        let pdf = &self.density[0];

        // Select a bin to sample from
        let sampled_bins = &self.per_bin_values[0];
        let sampled_bin = sampled_bins[pdf.sample(&mut rng)];

        // Generate a random data point using this PDF
        sampled_bin.sample(&mut rng)
    }
}

impl CustomShortErrorProfile {
    pub fn new(model_params: &encoding::ErrorModelParams, seed: Option<u64>) -> Self {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Start generating the quality score distributions
        //let mut scores = Vec::new();
        //let mut distributions = Vec::new();

        //println!(
        //    "Generating quality score distributions for {} positions",
        //    model_params.binned_quality_density.len()
        //);

        /*
        // Simulate quality score at each position in the read
        for i in 0..model_params.binned_quality_density.len() {
            let bins = &model_params.binned_quality_density[i];

            //let mut centers = Vec::new();
            //let mut counts = Vec::new();

            // Since we're not using individual data points, calculate centers for each bin
            for (bndx, count) in bins.iter().enumerate() {
                centers.push(if bndx == 0 {
                    model_params.bin_size as f64 / 2.0
                } else {
                    (((bndx - 1) * model_params.bin_size as usize)
                        + (bndx * model_params.bin_size as usize)) as f64
                        / 2.0
                });
                // We store center * count values, since that's the amount of observations in that bin
                for _ in 0..*count {
                    // We generate random values within this bin. We could also use bin centers
                    // but I think this generates more realistic scores
                    let val = if bndx == 0 {
                        rng.gen_range(0..model_params.bin_size) as u32
                    } else {
                        rng.gen_range(
                            ((bndx - 1) * model_params.bin_size as usize)
                                ..(bndx * model_params.bin_size as usize),
                        ) as u32
                    };
                    // aren't really centers but w/e
                    counts.push(val as f64);
                }
            }

            // Estimate the bandwidth parameter using Silverman's rule of thumb, not the best but it's fast
            let bandwidth =
                1.06 * stddev(counts.clone().into_iter()) * (counts.len() as f64).powf(-1.0 / 5.0);
            //let bandwidth = 1.06
            //    * stddev(centers.clone().into_iter())
            //    * (centers.len() as f64).powf(-1.0 / 5.0);

            distributions.push(QualityScoreDistribution { bandwidth, centers });
        }
            */

        Self {
            model_params: model_params.clone(),
            quality_score_pdf: CustomPDF::new(&model_params.binned_quality_density),
            read_length_pdf: CustomPDF::new(&vec![model_params.read_length_bins.clone()]),
            insert_size_pdf: if model_params.insert_size_bins.is_some() {
                Some(CustomPDF::new(&vec![model_params
                    .insert_size_bins
                    .as_ref()
                    .unwrap()
                    .clone()]))
            } else {
                None
            },
        }
    }
}

impl base::ErrorProfile for CustomShortErrorProfile {
    /**
     * Return a random read length based on the distribution of read lengths
     * specified by the model parameters.
     */
    fn get_read_length(&self, seed: Option<u64>) -> u16 {
        //let mut rng = match seed {
        //    Some(s) => StdRng::seed_from_u64(s),
        //    None => StdRng::from_entropy(),
        //};

        self.read_length_pdf.sample(seed) as u16

        // Sample a new value, truncate the resulting float into a u16
        //rng.sample(
        //    &Normal::new(
        //        self.model_params.read_length_mean,
        //        self.model_params.read_length_std,
        //    )
        //    .unwrap(),
        //)
        //.floor() as u16
    }

    /**
     * Return a random insert size based on the distribution of insert sizes
     * specified by the model parameters.
     */
    fn get_insert_size(&self, seed: Option<u64>) -> u16 {
        //let mut rng = match seed {
        //    Some(s) => StdRng::seed_from_u64(s),
        //    None => StdRng::from_entropy(),
        //};

        if self.insert_size_pdf.is_none() {
            0
        } else {
            self.insert_size_pdf.as_ref().unwrap().sample(seed) as u16
        }

        // Sample a new value, truncate the resulting float into a u16
        //rng.sample(
        //    &Normal::new(
        //        self.model_params.insert_size_mean,
        //        self.model_params.insert_size_std,
        //    )
        //    .unwrap(),
        //)
        //.floor() as u16
    }

    /**
     * Generates a random read length using a normal distribution.
     */
    fn get_random_read_length(&self, seed: Option<u64>) -> u16 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Sample a new value, truncate the resulting float into a u16
        rng.sample(
            &Normal::new(
                self.model_params.read_length_mean,
                self.model_params.read_length_std,
            )
            .unwrap(),
        )
        .floor() as u16
    }

    /**
     * Generates a random insert size using a normal distribution.
     */
    fn get_random_insert_size(&self, seed: Option<u64>) -> u16 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Sample a new value, truncate the resulting float into a u16
        rng.sample(
            &Normal::new(
                self.model_params.insert_size_mean,
                self.model_params.insert_size_std,
            )
            .unwrap(),
        )
        .floor() as u16
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
        //println!("Simulating phred scores for read of length {}", seq_length);
        //println!("qualities length: {}", self.model_params.qualities.len());

        let mut scores = Vec::new();

        // Generate a random quality score for each position in the read
        for i in 0..seq_length {
            // If we don't have any PDFs for this position, just use the last position PDF
            let seq_pos = if i >= self.quality_score_pdf.density.len() {
                self.quality_score_pdf.density.len() - 1
            } else {
                i
            };

            // Sample a random quality score from the PDF, truncates to a u8 but this should be
            // fine
            scores.push(self.quality_score_pdf.sample_with_index(seq_pos, seed) as u8);
        }

        return scores;

        /*
        for i in 0..seq_length {
            let distribution = if i < self.quality_score_distributions.len() {
                &self.quality_score_distributions[i]
            } else {
                &self.quality_score_distributions[self.quality_score_distributions.len() - 1]
            };

            let qs = distribution.centers.choose(&mut rng).unwrap();
            let normal = Normal::new(*qs, distribution.bandwidth).unwrap();
            let sim_qs = normal.sample(&mut rng);

            scores.push(sim_qs.round() as u8);
        }
        */
        //return scores;

        /*
        // Simulate quality score at each position in the read
        for i in 0..seq_length {
            println!("Simulating score for position {}", i);
            // Get bins for this position
            let bins = if i > self.model_params.qualities.len() {
                &self.model_params.qualities[self.model_params.qualities.len() - 1]
            } else {
                &self.model_params.qualities[i]
            };

            let mut centers = Vec::new();

            // Since we're not using individual data points, calculate centers for each bin
            for (bndx, count) in bins.iter().enumerate() {
                let sub_count = ((*count as f64) * 0.50) as usize;
                //println!("Bin {} has {} counts", bndx, count);
                // We store center * count values, since that's the amount of observations in that bin
                //for _ in 0..*count {
                for _ in 0..sub_count {
                    // We generate random values within this bin. We could also use bin centers
                    // but I think this generates more realistic scores
                    let val = if bndx == 0 {
                        rng.gen_range(0..self.model_params.bin_size) as u32
                    } else {
                        rng.gen_range(
                            ((bndx - 1) * self.model_params.bin_size as usize)
                                ..(bndx * self.model_params.bin_size as usize),
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
            //println!("scores {:?}", scores);
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
        */

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

        // Convert the vector of k-mer probabilities into a hashmap for faster lookups
        let kmer_probs = self
            .model_params
            .probabilities
            .iter()
            .map(|(k, p)| (k.clone(), p.clone()))
            .collect::<HashMap<u32, Vec<(u32, f32)>>>();
        //let mut kmer_probs = HashMap::new();

        //for (kmer, probs) in self.model_params.probabilities.iter() {
        //    kmer_probs.insert(kmer.clone(), probs.clone());
        //}

        for i in 0..sequence.len() {
            if (i + self.model_params.kmer_size) > sequence.len() {
                break;
            }

            // Grab the kmer
            let kmer = &new_sequence[i..i + self.model_params.kmer_size];
            // Encode it. If it can't be encoded b/c of characters other than ACGT, skip it
            let encoded_kmer =
                match encoding::three_bit_encode_kmer(kmer, self.model_params.kmer_size) {
                    Ok(k) => k,
                    Err(_) => continue,
                };

            // Find the list of alternate kmers, if this kmer hasn't been observed then move on
            //let alts = match self.model_params.probabilities.get(&encoded_kmer) {
            let alts = match kmer_probs.get(&encoded_kmer) {
                Some(vs) => vs,
                None => continue,
            };

            // Weighted distribution sampling
            let weighted_dist =
                WeightedAliasIndex::new(alts.iter().map(|(_, w)| *w).collect()).unwrap();

            // Get an alternate kmer which could include the og kmer itself
            let alt_encoded_kmer = alts[weighted_dist.sample(&mut rng)].0;

            // Decode the kmer into a byte array
            let alt_kmer = encoding::three_bit_decode_kmer(
                alt_encoded_kmer,
                self.model_params.kmer_size,
                true,
            )
            .unwrap();

            // Replace this kmer with the alt in the sequence. Splice is neat in that if alt_kmer has deletions,
            // it will delete the elements at those indices from new_sequence too.
            new_sequence.splice(i..i + self.model_params.kmer_size, alt_kmer);
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

    /**
     * Uses mean read lengths and insert size even though the minimum genome size necessary
     * could be larger than that due to RNG.
     */
    fn minimum_genome_size(&self) -> u16 {
        //2 * self.get_read_length() + self.get_insert_size()
        (2.0 * self.model_params.read_length_mean + self.model_params.insert_size_mean) as u16
    }

    fn is_long_read(&self) -> bool {
        self.model_params.is_long
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
    */
}
