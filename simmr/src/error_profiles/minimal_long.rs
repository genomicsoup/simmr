/**
 * file: minimal_long.rs
 * desc: The minimal long-read error profile. This generates phred scores using a normal
 * distribution, generates read lengths usning a gamma distribution, and assumes a uniform
 * probability of point mutations.
 *       
 */
use rand::prelude::SliceRandom;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Gamma, Normal};

use super::base;
use crate::util;

pub struct MinimalLongErrorProfile {}

impl base::ErrorProfile for MinimalLongErrorProfile {
    /**
     * Unusable functions. These can't and shouldn't be used for long-read error profiles.
     * Calling these when simulating long reads (which should never happen) will cause a panic.
     */
    fn get_read_length(&self) -> u16 {
        panic!("get_read_length() is not usable when simulating minimal long reads");
    }

    fn get_insert_size(&self) -> u16 {
        panic!("get_insert_size() is not usable when simulating minimal long reads");
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

        // Convert mean and standard deviation to shape and scale parameters
        let shape = (mean / std_dev).powf(2.0);
        let scale = std_dev.powf(2.0) / mean;

        // Sample a new value, truncate the resulting float into a u16
        rng.sample(&Gamma::new(shape, scale).unwrap()).floor() as u16
    }

    /**
     * Simulates phred scores using a normal distribution. The mean phred score is 30.
     */
    fn simulate_phred_scores(&self, seq_length: usize, seed: Option<u64>) -> Vec<u8> {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Assume mean accuracy 0.99 (phred score of 20) and a std. dev of 0.01, for long reads
        // which tend have worse quality scores this should probably be lowered
        let normal = Normal::new(util::convert_phred_to_accuracy(20), 0.05).unwrap();

        (0..seq_length)
            .into_iter()
            .map(|_| rng.sample(&normal))
            //.map(util::convert_probability_to_phred)
            // Prevent accuracies > 1.0
            .map(|acc| acc.min(0.9999))
            .map(util::convert_accuracy_to_phred)
            .collect::<Vec<u8>>()
    }

    /**
     * Simulates random point mutations using the generated quality scores. If a score is less than
     * a random probability, the nucleotide will get mutated into another one. Mutations have a uniform
     * chance of mutating into any other nucleotide.
     */
    fn simulate_point_mutations(
        &self,
        sequence: &[u8],
        quality: &Vec<u8>,
        seed: Option<u64>,
    ) -> Vec<u8> {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Sequence and quality vectors should be the same size, if not something is wrong.
        // Simulate point mutations.
        sequence
            .iter()
            .zip(quality.iter())
            .into_iter()
            .map(|(nt, q)| {
                // This is higher than the current phred score, so mutate
                if rng.gen::<f32>() > util::convert_phred_to_accuracy(*q) {
                    match nt {
                        b'A' => [b'C', b'G', b'T'].choose(&mut rng).unwrap(),
                        b'C' => [b'A', b'G', b'T'].choose(&mut rng).unwrap(),
                        b'T' => [b'A', b'C', b'G'].choose(&mut rng).unwrap(),
                        b'G' => [b'A', b'C', b'T'].choose(&mut rng).unwrap(),
                        // shouldn't ever get here
                        x => x,
                    }
                } else {
                    nt
                }
            })
            .map(|nt| *nt)
            .collect::<Vec<u8>>()
    }

    /**
     * This doesn't actually do anything for minimal profiles but since it's used by other long
     * read profiles, it still needs to return the sequence.
     */
    fn simulate_errors(&self, sequence: &[u8], seed: Option<u64>) -> Vec<u8> {
        sequence.to_vec()
    }

    // TODO: Hardcoded to the mean read length for long reads, should make this
    // user configurable
    fn minimum_genome_size(&self) -> u16 {
        20_000
    }

    fn is_long_read(&self) -> bool {
        true
    }
}
