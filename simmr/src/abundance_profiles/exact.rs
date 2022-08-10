/**
 * file: exact.rs
 * desc: The exact abundance profile generates exactly N reads per genome, it doesn't
 *       consider abundance values or account for size.
 */
use super::base;

use crate::genome;

pub struct ExactAbundanceProfile {}

impl base::AbundanceProfile for ExactAbundanceProfile {
    fn is_size_aware(&self) -> bool {
        false
    }

    fn determine_abundances(&self, total_reads: usize, num_genomes: usize) -> Vec<(usize, f64)> {
        let abundances = 100.0 / (num_genomes as f64);

        (0..num_genomes)
            .into_iter()
            .map(|_| (total_reads, abundances))
            .collect()
    }

    fn adjust_for_size(
        &self,
        _genomes: &Vec<genome::Genome>,
        read_abundances: &Vec<(usize, f64)>,
        _read_length: usize,
        _paired: bool,
    ) -> Vec<(usize, f64)> {
        // This abundance profile can't be adjusted for genome size
        read_abundances.clone()
    }
}
