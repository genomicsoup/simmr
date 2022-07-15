/**
 * file: simulate.rs
 * desc: Simulate reads.
 */
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::path;
use std::sync::atomic::{AtomicU32, Ordering};

use crate::abundance_profiles;
use crate::error_profiles;
use crate::genome;
use crate::util;

/**
 * STRUCTS
 */

/**
 * Models a single read which is a sequence of nucleotides and their corresponding quality scores.
 */
#[derive(Debug)]
pub struct SingleRead {
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

/**
 * Models a single PE read consisting of a unique ID and forward/reverse reads.
 */
#[derive(Debug)]
pub struct PERead {
    pub id: u32,
    pub forward: SingleRead,
    pub reverse: SingleRead,
}

/**
 * Models a single long read
 */
#[derive(Debug)]
pub struct LongRead {
    pub id: u32,
    pub read: SingleRead,
}

/**
 * FUNCTIONS
 */

/**
 * Keeps track of simulated read IDs across genomes. Returns and increments the current read ID
 * when a new simulated read is generated.
 */
fn get_read_id() -> u32 {
    static READ_ID: AtomicU32 = AtomicU32::new(0);

    READ_ID.fetch_add(1, Ordering::SeqCst)
}

/**
 * Simulate reads using the given genomes and error profile. This function also determines
 * genome abundances using the provided abundance profile.
 *
 * args
 *  num_reads:         total number of reads to generate
 *  genomes:           genomes to use
 *  error_profile:     error profile to use when simulating sequencing errors
 *  abundance_profile: genome abundance profile
 *  seed:              an optional seed for reproducibility
 *
 * returns
 *  a vector containing simulation metadata and the simulated reads, each element:
 *      0: path to the original genome
 *      1: generated unique ID for the genome
 *      2: number of reads simulated from this genome
 *      3: genome abundance
 *      4: simulated reads
 */
pub fn simulate_pe_reads<
    T: error_profiles::ErrorProfile + ?Sized,
    U: abundance_profiles::AbundanceProfile + ?Sized,
>(
    num_reads: usize,
    genomes: &Vec<genome::Genome>,
    error_profile: &T,
    abundance_profile: &U,
    seed: Option<u64>,
) -> Vec<(path::PathBuf, genome::UUID, usize, f64, Vec<PERead>)> {
    // Figure out abundance levels
    //let abundances = abundance_profile.determine_abundances(num_reads, &genomes);
    let abundances = abundance_profile.determine_abundances(num_reads, genomes.len());
    let mut all_reads = Vec::new();

    for ((genome_reads, abundance), genome) in abundances.iter().zip(genomes.iter()) {
        let reads =
            simulate_pe_reads_from_genome(*genome_reads, genome, error_profile, seed).unwrap();

        // Collect+asociate reads with genome metadata for output
        all_reads.push((
            genome.filepath.clone(),
            genome.uuid.clone(),
            *genome_reads,
            *abundance,
            reads,
        ));
    }

    all_reads
}

/**
 * Simulate a set of paired end reads, up to num_reads, using sequences from the given genome, and
 * according to an error profile.
 *
 * args
 *  num_reads:     the number of PE reads to generate
 *  genome:        genome to use for simulation
 *  error_profile: error profile to use when simulating sequencing errors
 *  seed:          an optional seed for reproducibility
 *
 * returns
 *  a result containing either an error string or the set of simulated PE reads
 */
pub fn simulate_pe_reads_from_genome<T: error_profiles::ErrorProfile + ?Sized>(
    num_reads: usize,
    genome: &genome::Genome,
    error_profile: &T,
    seed: Option<u64>,
) -> Result<Vec<PERead>, String> {
    let mut reads = Vec::new();
    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    for _ in 0..num_reads {
        // Select a random sequence from the genome
        //let i = rng.gen_range(0..genome.num_seqs);
        let seq = genome.sequence[rng.gen_range(0..genome.num_seqs)].clone();
        //let seq = genome.sequence[i].clone();
        let pe_seed: u64 = rng.gen();

        reads.push(simulate_pe_read(&seq, error_profile, Some(pe_seed)).unwrap())
    }

    Ok(reads)
}

/**
 * Simulate a single paired end read using the given sequence from a single genome and according
 * to an error profile and distribution.
 *
 * args
 *  sequence:      a sequence record from a genome
 *  error_profile: error profile to use when simulating sequencing errors
 *  read_id:       generated and unique read ID to assign to the simulated read
 *  seed:          an optional seed for reproducibility
 *
 * returns
 *  a result containing either an error or the simulated PE read
 */
pub fn simulate_pe_read<T: error_profiles::ErrorProfile + ?Sized>(
    //genome: &genome::Genome,
    sequence: &genome::Seq,
    error_profile: &T,
    seed: Option<u64>,
) -> Result<PERead, String> {
    let read_length = error_profile.get_read_length() as usize;
    let insert_size = error_profile.get_insert_size() as usize;

    // Genome size must be at min. this length to generate a PE read
    let required_size = 2 * read_length + insert_size;

    // We should never get this error b/c prior to calling this function we remove sequences
    // that are smaller than needed but this is here just in case
    if sequence.size <= required_size {
        return Err(format!(
            "Genome size ({}nt) is smaller than the required length ({})",
            sequence.size, required_size
        ));
    }

    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    // Forward read starts and ends
    let fwd_start: usize = rng.gen_range(0..(sequence.size - required_size));
    let fwd_end = fwd_start + read_length;

    // If these are out of bounds (first is reverse read end and second is it's start position),
    // regenerate the start and end positions for the reverse read
    let (rev_end, rev_start) = if (fwd_end + insert_size) >= sequence.size
        || (fwd_end + insert_size + read_length) >= sequence.size
    {
        // TODO: This has the side effect of not respecting the insert size, fix this at some point
        let new_rev_end = rng.gen_range(fwd_start..(sequence.size - required_size));

        (new_rev_end, new_rev_end + read_length)
    } else {
        (fwd_end + insert_size, fwd_end + insert_size + read_length)
    };

    // Generate the sequences from the genome
    let mut fwd_read = sequence.seq[fwd_start..fwd_end].to_owned();
    let mut rev_read = sequence.seq[rev_end..rev_start].to_owned();

    // Generate quality scores for both reads
    let fwd_quality = error_profile.simulate_phred_scores(read_length, seed);
    let rev_quality = error_profile.simulate_phred_scores(read_length, seed);

    // Simulate point mutations using this profile's error distribution
    fwd_read = error_profile.simulate_point_mutations(&fwd_read, &fwd_quality, seed);
    rev_read = error_profile.simulate_point_mutations(&rev_read, &rev_quality, seed);

    // Finally, fill out the read struct
    let pe_read = PERead {
        id: get_read_id(),
        forward: SingleRead {
            sequence: fwd_read,
            quality: fwd_quality,
        },
        reverse: SingleRead {
            // Don't forget the reverse read must be the reverse complement since it is sequenced
            // from 3' -> 5'
            // TODO: move this up before mutation occurs?
            sequence: util::reverse_complement(&rev_read),
            quality: rev_quality,
        },
    };

    Ok(pe_read)
}

/**
 * Simulate long reads using the given genomes and error profile. This function also determines
 * genome abundances using the provided abundance profile.
 *
 * args
 *  num_reads:         total number of reads to generate
 *  genomes:           genomes to use
 *  error_profile:     error profile to use when simulating sequencing errors
 *  abundance_profile: genome abundance profile
 *  seed:              an optional seed for reproducibility
 *
 * returns
 *  a vector containing simulation metadata and the simulated reads, each element:
 *      0: path to the original genome
 *      1: generated unique ID for the genome
 *      2: number of reads simulated from this genome
 *      3: genome abundance
 *      4: simulated reads
 */
pub fn simulate_long_reads<
    T: error_profiles::ErrorProfile + ?Sized,
    U: abundance_profiles::AbundanceProfile + ?Sized,
>(
    num_reads: usize,
    genomes: &Vec<genome::Genome>,
    error_profile: &T,
    abundance_profile: &U,
    seed: Option<u64>,
) -> Vec<(path::PathBuf, genome::UUID, usize, f64, Vec<LongRead>)> {
    // Figure out abundance levels
    let abundances = abundance_profile.determine_abundances(num_reads, genomes.len());
    let mut all_reads = Vec::new();

    for ((genome_reads, abundance), genome) in abundances.iter().zip(genomes.iter()) {
        let reads =
            simulate_long_reads_from_genome(*genome_reads, genome, error_profile, seed).unwrap();

        // Collect+asociate reads with genome metadata for output
        all_reads.push((
            genome.filepath.clone(),
            genome.uuid.clone(),
            *genome_reads,
            *abundance,
            reads,
        ));
    }

    all_reads
}

/**
 * Simulate a set of long reads, up to num_reads, using sequences from the given genome, and
 * according to an error profile.
 *
 * args
 *  num_reads:     the number of PE reads to generate
 *  genome:        genome to use for simulation
 *  error_profile: error profile to use when simulating sequencing errors
 *  seed:          an optional seed for reproducibility
 *
 * returns
 *  a result containing either an error string or the set of simulated PE reads
 */
pub fn simulate_long_reads_from_genome<T: error_profiles::ErrorProfile + ?Sized>(
    num_reads: usize,
    genome: &genome::Genome,
    error_profile: &T,
    seed: Option<u64>,
) -> Result<Vec<LongRead>, String> {
    let mut reads = Vec::new();
    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };
    let mut simulated_read_count = 0;

    while simulated_read_count < num_reads {
        // Select a random sequence from the genome
        let seq = genome.sequence[rng.gen_range(0..genome.num_seqs)].clone();
        let read_seed: u64 = rng.gen();
        // TODO: Will have to do something about this b/c this could cause an infinite loop
        // under certain edge cases
        let simulated = match simulate_long_read(&seq, error_profile, Some(read_seed)) {
            Ok(r) => r,
            Err(e) => continue,
        };

        reads.push(simulated);
    }

    Ok(reads)
}

/**
 * Simulate a single long read using the given sequence from a single genome and according
 * to an error profile and distribution.
 *
 * args
 *  sequence:      a sequence record from a genome
 *  error_profile: error profile to use when simulating sequencing errors
 *  seed:          an optional seed for reproducibility
 *
 * returns
 *  a result containing either an error or the simulated read
 */
pub fn simulate_long_read<T: error_profiles::ErrorProfile + ?Sized>(
    sequence: &genome::Seq,
    error_profile: &T,
    seed: Option<u64>,
) -> Result<LongRead, String> {
    // The genome might not be long enough for our random read length, so try a few times and
    // if doesn't work out move to another genome.
    let mut read_length = error_profile.get_random_read_length(seed);

    if sequence.size <= read_length as usize {
        for _ in 0..10 {
            read_length = error_profile.get_random_read_length(seed);

            if read_length as usize > sequence.size {
                break;
            }
        }
    }

    // It's possible the size is still smaller than the read length, if so we'll return an
    // error and try a new genome
    if sequence.size <= read_length as usize {
        return Err(format!(
            "Genome size ({}nt) is smaller than the read length ({})",
            sequence.size, read_length
        ));
    }

    let mut rng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    // Start and end positions of the read
    let read_start: usize = rng.gen_range(0..read_length as usize);
    let mut read_end: usize = read_start + read_length as usize;

    // if the end of the read is out of bounds, then regenerate the end position
    if read_end >= sequence.size {
        read_end = rng.gen_range(read_start..sequence.size);
        read_length = (read_end - read_start) as u16;
    }

    // Generate the sequence from the genome
    let mut read_sequence = sequence.seq[read_start..read_end].to_owned();

    // Generate quality scores
    let read_quality = error_profile.simulate_phred_scores(read_length as usize, seed);

    // Simulate errors, this only modifies the sequence if using a custom profile
    read_sequence = error_profile.simulate_errors(&read_sequence, seed);

    // Simulate point mutations, this only modifies the sequence if using a minimal profile
    read_sequence = error_profile.simulate_point_mutations(&read_sequence, &read_quality, seed);

    // Finally, fill out the read struct
    let long_read = LongRead {
        id: get_read_id(),
        read: SingleRead {
            sequence: read_sequence,
            quality: read_quality,
        },
    };

    Ok(long_read)
}

#[cfg(test)]
#[path = "tests/simulate_tests.rs"]
mod simulate_tests;
