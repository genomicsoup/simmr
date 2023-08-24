/**
 * file: custom.rs
 * desc: A custom abundance profile which allows the user to specify the abundance of each
 *       genome.
 */
use super::base;

use crate::genome;

pub struct CustomAbundanceProfile {
    pub size_adjusted: bool,
    pub abundances: Vec<f64>,
}

impl base::AbundanceProfile for CustomAbundanceProfile {
    fn is_size_aware(&self) -> bool {
        self.size_adjusted
    }

    fn determine_abundances(
        &self,
        total_reads: usize,
        num_genomes: usize,
        //genomes: &Vec<genome::Genome>,
    ) -> Vec<(usize, f64)> {
        let total_abundance: f64 = self.abundances.iter().sum();
        // abundances should sum up to roughly 1.0, we give some leeway. if it doesn't sum, then we
        // have to normalize it
        if total_abundance < 0.99 || total_abundance > 1.01 {
            self.abundances
                .iter()
                .map(|abund| {
                    (
                        (total_reads as f64 * (abund / total_abundance)).ceil() as usize,
                        abund / total_abundance,
                    )
                })
                .collect()
        } else {
            self.abundances
                .iter()
                .map(|abund| ((total_reads as f64 * abund).ceil() as usize, *abund))
                .collect::<Vec<_>>()
        }
    }

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
        //genomes
        //    .into_iter()
        //    .zip(read_abundances)
        //    .map(|(g, (nreads, abund))| {
        //        (
        //            ((total_coverage * total_reads * (*abund / 100.0))
        //                / base::coverage(*nreads, read_length, g.size, true))
        //            .ceil() as usize,
        //            *abund,
        //        )
        //    })
        //    .collect()

        let total_adjusts: f64 = genomes
            .into_iter()
            .zip(read_abundances)
            .map(|(g, (_, abund))| g.size as f64 * *abund)
            .sum();

        genomes
            .into_iter()
            .zip(read_abundances)
            .map(|(g, (nreads, abund))| {
                (
                    (total_reads * ((*abund * g.size as f64) / total_adjusts)).ceil() as usize,
                    *abund,
                )
            })
            .collect()
    }
}
