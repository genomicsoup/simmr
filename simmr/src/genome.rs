/**
 * file: genome.rs
 * desc: Model a genome and its sequence.
 */
use needletail::parse_fastx_file;
use needletail::Sequence;
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
    pub uuid: UUID,
    pub filepath: path::PathBuf,
    pub sequence: Vec<Seq>,
    pub size: usize,
    pub num_seqs: usize,
    pub abundance: Option<f64>,
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
    pub fn from_fasta(filepath: &str) -> Result<Genome, String> {
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

        let genome = Genome {
            // this is broken so just disable it for now and replace it with something gross
            uuid: UUID::from(util::generate_id(None)),
            filepath: path::PathBuf::from(filepath),
            sequence: sequences.clone(),
            size: sequences.iter().map(|s| s.seq.len()).sum(),
            num_seqs: sequences.len(),
            abundance: None,
        };

        Ok(genome)
    }
}

#[cfg(test)]
#[path = "tests/genome_tests.rs"]
mod genome_tests;
