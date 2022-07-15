/**
 * file: uniform.rs
 * desc: The uniform abundance profile generates uniform abundances across all genomes.
 *       This is based on the number of reads to generate.
 */
use super::base;

pub struct UniformAbundanceProfile {}

impl base::AbundanceProfile for UniformAbundanceProfile {
    fn determine_abundances(
        &self,
        num_reads: usize,
        num_genomes: usize, //genomes: &Vec<genome::Genome>,
    ) -> Vec<(usize, f64)> {
        // Uniform distribution splits abundances evenly across genomes
        let reads_per_genome = ((num_reads as f64) / (num_genomes as f64)).ceil() as usize;
        let abundances = 100.0 / (num_genomes as f64);

        (0..num_genomes)
            .into_iter()
            .map(|_| (reads_per_genome, abundances))
            .collect()
    }
}
