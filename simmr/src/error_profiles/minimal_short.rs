/**
 * file: minimal_short.rs
 * desc: The minimal error profile generates phred scores using a normal distribution and
 *       assumes a uniform probability of point mutations. Insert sizes and read lengths are
 *       also generated using a normal distribution.
 */
use rand::prelude::SliceRandom;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::Normal;

use super::base;
use crate::util;

pub struct MinimalShortErrorProfile {
    pub read_length: u16,
    pub insert_size: u16,
    pub mean_phred_score: u8,
    pub insert_size_std: f64,
    pub read_length_std: f64,
}

impl base::ErrorProfile for MinimalShortErrorProfile {
    /**
     * Unusable functions. These can't and shouldn't be used for short read error profiles.
     * Calling these when simulating short reads (which should never happen) will cause a panic.
     */
    fn simulate_errors(&self, sequence: &[u8], seed: Option<u64>) -> Vec<u8> {
        panic!("simulate_errors() is not usable when simulating basic short reads");
    }

    fn get_read_length(&self, seed: Option<u64>) -> u16 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Sample a new value, truncate the resulting float into a u16
        rng.sample(&Normal::new(self.read_length as f64, self.read_length_std).unwrap())
            .floor() as u16
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
        rng.sample(&Normal::new(self.read_length as f64, self.read_length_std).unwrap())
            .floor() as u16
    }

    fn get_insert_size(&self, seed: Option<u64>) -> u16 {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Sample a new value, truncate the resulting float into a u16
        rng.sample(&Normal::new(self.insert_size as f64, self.insert_size_std).unwrap())
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
        rng.sample(&Normal::new(self.insert_size as f64, self.insert_size_std).unwrap())
            .floor() as u16
    }

    fn simulate_phred_scores(&self, seq_length: usize, seed: Option<u64>) -> Vec<u8> {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        // Assume mean quality score of 30 and std. dev of 10
        let normal = Normal::new(self.mean_phred_score as f32, 10.0).unwrap();
        //let normal = Normal::new(util::convert_phred_to_accuracy(15), 0.01).unwrap();

        (0..seq_length)
            .into_iter()
            .map(|_| rng.sample(&normal))
            //.map(util::convert_probability_to_phred)
            // Prevent accuracies > 1.0
            //.map(|acc| acc.min(0.9999))
            //.map(util::convert_accuracy_to_phred)
            .map(|p| p.floor() as u8)
            .collect::<Vec<u8>>()
    }

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

        // These should be the same size, if not something is wrong. Simulate point mutations
        // using a uniform distribution
        let s = sequence
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
            .collect::<Vec<u8>>();

        s
    }

    fn minimum_genome_size(&self) -> u16 {
        //2 * self.get_read_length() + self.get_insert_size()
        2 * self.read_length + self.insert_size
    }

    fn is_long_read(&self) -> bool {
        false
    }
}
