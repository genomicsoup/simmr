/**
 * file: base.rs
 * desc: Base abundance profile trait which all abundance models inherit from.
 */
pub trait AbundanceProfile {
    // Figure out abundance levels for each genome
    fn determine_abundances(
        &self,
        num_reads: usize,
        num_genomes: usize, //genomes: &Vec<genome::Genome>,
    ) -> Vec<(usize, f64)>;
}

impl<T: ?Sized> AbundanceProfile for Box<T>
where
    T: AbundanceProfile,
{
    fn determine_abundances(
        &self,
        num_reads: usize,
        num_genomes: usize, //genomes: &Vec<genome::Genome>,
    ) -> Vec<(usize, f64)> {
        (**self).determine_abundances(num_reads, num_genomes)
    }
}
