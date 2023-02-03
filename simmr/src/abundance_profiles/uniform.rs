/**
 * file: uniform.rs
 * desc: The uniform abundance profile generates uniform abundances across all genomes.
 */
use super::base;

use crate::{abundance_profiles::base::coverage, genome};

pub struct UniformAbundanceProfile {
    pub size_aware: bool,
}

impl base::AbundanceProfile for UniformAbundanceProfile {
    fn is_size_aware(&self) -> bool {
        self.size_aware
    }

    fn determine_abundances(
        &self,
        total_reads: usize,
        num_genomes: usize,
        //genomes: &Vec<genome::Genome>,
    ) -> Vec<(usize, f64)> {
        //let num_genomes = genomes.len();
        // First normalize based on genome sizes
        //let total_size =
        // Uniform distribution splits abundances evenly across genomes
        let reads_per_genome = ((total_reads as f64) / (num_genomes as f64)).ceil() as usize;
        let abundances = 100.0 / (num_genomes as f64);

        (0..num_genomes)
            .into_iter()
            .map(|_| (reads_per_genome, abundances))
            .collect()
    }

    //fn adjust_for_size(
    //    &self,
    //    total_coverage: f64,
    //    genome_coverage: f64,
    //    total_reads: usize,
    //    abundance: f64,
    //) -> f64 {
    //    (**self).adjust_for_size(total_coverage, genome_coverage, total_reads, abundance)
    //}
    fn adjust_for_size(
        &self,
        genomes: &Vec<genome::Genome>,
        read_abundances: &Vec<(usize, f64)>,
        read_length: usize,
        paired: bool,
    ) -> Vec<(usize, f64)> {
        // Determine the sum of coverages for all genomes
        let total_coverage = base::total_coverage(
            genomes,
            &read_abundances.iter().map(|(n, _)| *n).collect(),
            read_length,
            paired,
        );
        let total_reads: f64 = read_abundances.iter().map(|(n, _)| *n as f64).sum();

        // Adjust the number of reads to simulate based on genome size adjustments. The rel.
        // abundance will stay the same but reads per genome may change. See trait docstring
        // for more info on this calculation. Formula is
        // genome reads = total coverage * total reads * rel. abundance / genome coverage
        genomes
            .into_iter()
            .zip(read_abundances)
            .map(|(g, (nreads, abund))| {
                (
                    ((total_coverage * total_reads * (*abund / 100.0))
                        / base::coverage(*nreads, read_length, g.size, true))
                    .ceil() as usize,
                    *abund,
                )
            })
            .collect()
    }
}
