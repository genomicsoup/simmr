/**
 * file: base.rs
 * desc: Base error profile trait which all error models inherit from.
 */

pub trait ErrorProfile {
    // Returns the length of the simulated PE read -- this is only used for PE reads
    fn get_read_length(&self, seed: Option<u64>) -> u16;
    // Returns a random read length for a simulated read -- this is primarily used for long reads
    fn get_random_read_length(&self, seed: Option<u64>) -> u16;
    // Returns the length of the simulated insert size
    fn get_insert_size(&self, seed: Option<u64>) -> u16;
    // Returns a random insert size for a simulated read
    fn get_random_insert_size(&self, seed: Option<u64>) -> u16;
    // Simulate phred scores based on the error model
    fn simulate_phred_scores(&self, seq_length: usize, seed: Option<u64>) -> Vec<u8>;
    // Simulate point mutations in the given sequence using profile-specific distributions
    fn simulate_point_mutations(
        &self,
        sequence: &[u8],
        quality: &Vec<u8>,
        seed: Option<u64>,
    ) -> Vec<u8>;
    // This introduces substitutions and indels using a kmer based method. This is intended only
    // for custom error profiles
    fn simulate_errors(&self, sequence: &[u8], seed: Option<u64>) -> Vec<u8>;
    // Calculate the minimum genome size necessary to simulate reads using the given read length
    // and insert sizes
    fn minimum_genome_size(&self) -> u16;
    // Returns true if the error profile is used for generating long reads
    fn is_long_read(&self) -> bool;
}

impl<T: ?Sized> ErrorProfile for Box<T>
where
    T: ErrorProfile,
{
    fn get_read_length(&self, seed: Option<u64>) -> u16 {
        (**self).get_read_length(seed)
    }

    fn get_random_read_length(&self, seed: Option<u64>) -> u16 {
        (**self).get_random_read_length(seed)
    }

    fn get_insert_size(&self, seed: Option<u64>) -> u16 {
        (**self).get_insert_size(seed)
    }

    fn get_random_insert_size(&self, seed: Option<u64>) -> u16 {
        (**self).get_random_insert_size(seed)
    }

    fn simulate_phred_scores(&self, seq_length: usize, seed: Option<u64>) -> Vec<u8> {
        (**self).simulate_phred_scores(seq_length, seed)
    }

    fn simulate_point_mutations(
        &self,
        sequence: &[u8],
        quality: &Vec<u8>,
        seed: Option<u64>,
    ) -> Vec<u8> {
        (**self).simulate_point_mutations(sequence, quality, seed)
    }

    fn simulate_errors(&self, sequence: &[u8], seed: Option<u64>) -> Vec<u8> {
        (**self).simulate_errors(sequence, seed)
    }

    fn minimum_genome_size(&self) -> u16 {
        (**self).minimum_genome_size()
    }

    fn is_long_read(&self) -> bool {
        (**self).is_long_read()
    }
}
