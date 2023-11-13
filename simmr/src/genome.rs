/**
 * file: genome.rs
 * desc: Model a genome and its sequence.
 */
use needletail::parse_fastx_file;
use needletail::Sequence;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::path;

use crate::util;

#[derive(Debug, Clone)]
pub struct UUID(Vec<u8>);

#[derive(Debug, Clone)]
pub struct Seq {
    pub id: Vec<u8>,
    pub uuid: u64,
    pub seq: Vec<u8>,
    pub size: usize,
}

#[derive(Debug, Clone)]
pub struct Genome {
    // Unique identifier for this genome
    pub uuid: UUID,
    // Genome (FASTA) filepath
    pub filepath: path::PathBuf,
    // Sequence records
    pub sequence: Vec<Seq>,
    // Total genome size in bp
    pub size: usize,
    // Number of sequences (separat records) in this genome
    pub num_seqs: usize,
    // Abundance of this genome in simulated samples
    pub abundance: Option<f64>,
    // If true, don't treat sequences as separate records but as a single contiguous sequence
    pub contiguous: bool,
}

// Display impls to format for print and get a free to_string() impl
impl std::fmt::Display for UUID {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", std::str::from_utf8(&self.0).unwrap())
    }
}

// Equality implementation so we can compare UUIDs
impl PartialEq for UUID {
    fn eq(&self, other: &Self) -> bool {
        if self.0.len() != other.0.len() {
            return false;
        }

        self.0.iter().zip(other.0.iter()).all(|(a, b)| *a == *b)
    }
}

// Convenience methods for converting some types into a UUID
impl From<u64> for UUID {
    fn from(u: u64) -> Self {
        UUID(format!("{:x}", u).as_bytes().to_vec())
    }
}

impl From<String> for UUID {
    fn from(s: String) -> Self {
        UUID(s.as_bytes().to_vec())
    }
}

impl std::fmt::Display for Genome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "<Genome ({}) {}>",
            self.filepath.file_name().unwrap().to_str().unwrap(),
            self.uuid.to_string()
        )
    }
}

impl Genome {
    /**
     * Constructs a new Genome object from the given FASTA file.
     */
    pub fn from_fasta(filepath: &str, contiguous: bool) -> Result<Genome, String> {
        // Holds all records for this fasta
        let mut sequences: Vec<Seq> = Vec::new();
        // Record iterator, any parse errors get returned to the caller
        let mut fasta_reader = match parse_fastx_file(filepath) {
            Ok(r) => r,
            Err(e) => return Err(format!("{}", e)),
        };

        while let Some(record_wrap) = fasta_reader.next() {
            // Fail on any parse errors even if all other records are fine
            let record = match record_wrap {
                Ok(r) => r,
                Err(e) => return Err(format!("{}", e)),
            };
            // Fasta header
            let header = record.id();
            // Normalized sequences, removes softmasking, etc.
            let seq = record.normalize(false);

            sequences.push(Seq {
                id: header.to_vec(),
                uuid: util::generate_id(None),
                seq: seq.to_vec(),
                size: seq.to_vec().len(),
            })
        }

        let genome = if contiguous {
            Genome {
                // this is broken so just disable it for now and replace it with something gross
                uuid: UUID::from(util::generate_id(None)),
                filepath: path::PathBuf::from(filepath),
                sequence: vec![Seq {
                    id: b"whole genome".to_vec(),
                    uuid: util::generate_id(None),
                    seq: sequences
                        .iter()
                        .map(|s| s.seq.clone())
                        .flat_map(|s| s.into_iter().chain(std::iter::once(b'N')))
                        .collect::<Vec<u8>>(),
                    size: sequences.iter().map(|s| s.seq.len()).sum(),
                }],
                size: sequences.iter().map(|s| s.seq.len()).sum(),
                num_seqs: 1,
                abundance: None,
                contiguous,
            }
        } else {
            Genome {
                // this is broken so just disable it for now and replace it with something gross
                uuid: UUID::from(util::generate_id(None)),
                filepath: path::PathBuf::from(filepath),
                sequence: sequences.clone(),
                size: sequences.iter().map(|s| s.seq.len()).sum(),
                num_seqs: sequences.len(),
                abundance: None,
                contiguous,
            }
        };

        //let genome = Genome {
        //    // this is broken so just disable it for now and replace it with something gross
        //    uuid: UUID::from(util::generate_id(None)),
        //    filepath: path::PathBuf::from(filepath),
        //    sequence: sequences.clone(),
        //    size: sequences.iter().map(|s| s.seq.len()).sum(),
        //    num_seqs: sequences.len(),
        //    abundance: None,
        //    contiguous,
        //};

        Ok(genome)
    }

    /*
    Return a random sequence from this genome.
    If the genome is contiguous, this will return the entire genome as a single sequence.
    pub fn random_sequence(&self, seed: Option<u64>) -> &Seq {
        let mut rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };
        if !self.contiguous {
            let idx = rng.gen_range(0..self.sequence.len());
            &self.sequence[idx]
        } else {
            let seqs = &self
                .sequence
                .iter()
                .map(|s| s.seq)
                .flat_map(|s| s.iter().cloned().chain(std::iter::once(b'N')))
                .collect::<Vec<u8>>();

            &Seq {
                id: b"whole genome".to_vec(),
                uuid: util::generate_id(None),
                seq: seqs.clone(),
                size: seqs.len(),
            }
        }
    }
    */
}

#[cfg(test)]
#[path = "tests/genome_tests.rs"]
mod genome_tests;
