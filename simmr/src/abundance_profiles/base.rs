/**
 * file: base.rs
 * desc: Base abundance profile trait which all abundance models inherit from.
 */
use crate::genome;

/**
 * TRAITS
 */
pub trait AbundanceProfile {
    /**
    Do abundance calculations account for genome size?
    */
    fn is_size_aware(&self) -> bool;

    /**
     Determine the number of reads and relative abundance for each genome given a total number
     of reads.

     args
        total_reads: total number of reads being simulated
        num_genomes: number of genomes being used to simulate reads

    returns
        a list of tuples where each element has the following members
            0: the number of reads to generate for genome i
            1: the estimated relative abundance of genome i
     */
    fn determine_abundances(
        &self,
        total_reads: usize,
        //genomes: &Vec<genome::Genome>,
        num_genomes: usize,
    ) -> Vec<(usize, f64)>;

    // Adjust the number of simulated reads for the given relative abundance based on genome size
    // tc*tr*ra / cov
    //fn adjust_for_size(
    //    &self,
    //    total_coverage: f64,
    //    genome_coverage: f64,
    //    total_reads: usize,
    //    abundance: f64,
    //) -> f64;

    /**
    Adjusts the estimated read abundances to account for genome size. This uses a very simple
    relative abundance calculation based on genome coverage:
       rel. abundance = genome coverage / total coverage * (genome reads / total reads)
    If we rearrange the equation to solve for the number of genome reads we need to simulate to
    account for its size, the equation becomes:
       genome reads = total coverage * total reads * rel. abundance / genome coverage

    args
        genomes:         a list of genomes used for read simulation
        read_abundances: a read abundances to simulate per genome (num reads and abundance values)
        read_length:     Length of a single read

    returns
        the size-adjusted read abundance values
    */
    fn adjust_for_size(
        &self,
        genomes: &Vec<genome::Genome>,
        read_abundances: &Vec<(usize, f64)>,
        read_length: usize,
        paired: bool,
    ) -> Vec<(usize, f64)>;
}

/**
 * IMPLEMENTATIONS
 */

impl<T: ?Sized> AbundanceProfile for Box<T>
where
    T: AbundanceProfile,
{
    fn is_size_aware(&self) -> bool {
        (**self).is_size_aware()
    }

    fn determine_abundances(
        &self,
        total_reads: usize,
        //genomes: &Vec<genome::Genome>,
        num_genomes: usize,
    ) -> Vec<(usize, f64)> {
        (**self).determine_abundances(total_reads, num_genomes)
        //(**self).determine_abundances(num_reads, genomes)
    }

    //fn adjust_for_size(
    //    &self,
    //    total_coverage: f64,
    //    genome_coverage: f64,
    //    total_reads: usize,
    //    abundance: f64,
    //) -> usize {
    //    (**self).adjust_for_size(total_coverage, genome_coverage, total_reads, abundance)
    //}
    fn adjust_for_size(
        &self,
        genomes: &Vec<genome::Genome>,
        read_abundances: &Vec<(usize, f64)>,
        read_length: usize,
        paired: bool, //total_coverage: f64,
                      //genome_coverage: f64,
                      //total_reads: usize,
                      //abundance: f64,
    ) -> Vec<(usize, f64)> {
        (**self).adjust_for_size(genomes, read_abundances, read_length, paired)
    }
}

/**
 * FUNCTIONS
 */

/**
 Calculate coverage for the given genome.

 args
    - num_reads:   number of reads generated for this genome
    - read_length: length of an individual read
    - genome_size: size of the genome

 returns
    genome coverage
*/
pub fn coverage(num_reads: usize, read_length: usize, genome_size: usize, paired: bool) -> f64 {
    if paired {
        (num_reads as f64 * read_length as f64 * 2.0) / genome_size as f64
    } else {
        (num_reads as f64 * read_length as f64) / genome_size as f64
    }
}

/**
 Calculate total coverage for the given set of genomes.

 args
    - genomes:     list of genome objects to estimate coverage for
    - num_reads:   a list of the number of reads simulated for each genome, each element
                   corresponds to an element in the genome array.
    - read_length: length of an individual read

 returns
    total coverage aka the sum of all genome coverages
*/
pub fn total_coverage(
    genomes: &Vec<genome::Genome>,
    num_reads: &Vec<usize>,
    read_length: usize,
    paired: bool,
) -> f64 {
    genomes
        .iter()
        .zip(num_reads)
        .map(|(g, r)| coverage(*r, read_length, g.size, paired))
        .sum()
}
