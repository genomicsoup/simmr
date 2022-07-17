/**
 * file: perfect_short.rs
 * desc: The perfect error profile generates perfect short reads. No substitutions/indels and
 *       high, uniform phred scores.
 */
use super::base;

pub struct PerfectShortErrorProfile {
    pub read_length: u16,
    pub insert_size: u16,
}

impl base::ErrorProfile for PerfectShortErrorProfile {
    /**
     * Unusable functions. These can't and shouldn't be used for short read error profiles.
     * Calling these when simulating short reads (which should never happen) will cause a panic.
     */
    fn simulate_errors(&self, _sequence: &[u8], _seed: Option<u64>) -> Vec<u8> {
        panic!("simulate_errors() is not usable when simulating perfect short reads");
    }

    fn get_read_length(&self) -> u16 {
        self.read_length
    }

    fn get_random_read_length(&self, seed: Option<u64>) -> u16 {
        self.read_length
    }

    fn get_insert_size(&self) -> u16 {
        self.insert_size
    }

    fn simulate_phred_scores(&self, seq_length: usize, _seed: Option<u64>) -> Vec<u8> {
        std::iter::repeat(60).take(seq_length).collect()
    }

    fn simulate_point_mutations(
        &self,
        sequence: &[u8],
        _quality: &Vec<u8>,
        _seed: Option<u64>,
    ) -> Vec<u8> {
        // Since this error profile is used to generate perfect reads we don't mutate anything
        sequence.iter().cloned().collect::<Vec<u8>>()
    }

    fn is_long_read(&self) -> bool {
        false
    }
}
