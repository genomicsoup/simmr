/**
 * file: alignment.rs
 * desc: Functions related to parsing alignment data and SAM/BAM formats.
 */
use bincode;
use hex;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

use shared::encoding;
use shared::util;

/**
 * STRUCTS
 */

/**
 * Models metadata from an alignment record present in a SAM/BAM file.
 *
 * fields
 *  cigar:    CIGAR string, which keeps info about alignment (mis)matches, indels, etc.
 *  sequence: sequence of the aligned read
 *  md_tag:   encodes (mis)matches and deletions in the reference sequence
 */
#[derive(Deserialize, Serialize, PartialEq, Debug)]
pub struct AlignmentRecord {
    pub cigar: Vec<u8>,
    pub sequence: Vec<u8>,
    pub md_tag: Vec<u8>,
}

// Maps a bit-encoded kmer to its alternate kmers and counts of those kmers.
pub type EncodedKmerCountMap = HashMap<u32, HashMap<u32, u32>>;

/**
 * FUNCTIONS
 */

/**
 * Fully expand the given cigar string. e.g., 2M1I3M2D into MMIMMMDD
 *
 * args
 *  cigar: the original, unexpanded CIGAR string
 *
 * returns
 *  an expanded CIGAR string
 */
pub fn expand_cigar(cigar: &[u8]) -> Vec<u8> {
    let mut new_cigar = Vec::new();
    let mut cigar_iter = cigar.iter();

    while let Some(c) = cigar_iter.next() {
        // The run number for this set of matches, mismatches, indels, etc.
        let mut num = vec![*c];
        let mut next_char = cigar_iter.next();

        loop {
            match next_char {
                None => break,
                Some(nc) => {
                    // If it's an individual number, it's still part of the run number
                    if *nc >= b'0' && *nc <= b'9' {
                        num.push(*nc);
                    } else {
                        // Otherwise this is a letter and represents the run type
                        break;
                    }
                    next_char = cigar_iter.next();
                }
            }
        }

        // Convert the number array into the actual run number
        let runs = match std::str::from_utf8(&num).unwrap().parse::<u32>() {
            Err(e) => {
                panic!(
                    "CIGAR string ({}) is probably malformed: {}",
                    util::bytes_to_string(&cigar),
                    e
                );
            }
            Ok(i) => i,
        };

        let run_type = next_char.unwrap();

        // The next_char should contain the type of operation (match, indel, etc.), expand that
        for _ in 0..runs {
            new_cigar.push(*run_type);
        }
    }

    new_cigar
}

/**
 * Fully expand the given MD tag.
 * e.g., 2C1T^GC into [(M, M), (M, M), (N, C), (M, M), (D, G), (D, C)]
 * The expanded tag contains both the error type
 *      M = match, N = mismatch, D = deletion
 * and the nucleotide that replaced it.
 *
 * args
 *  md: the original, raw MD tag
 *
 * returns
 *  an expanded MD tag
 */
pub fn expand_md_tag(md: &[u8]) -> Vec<(u8, u8)> {
    let mut new_md: Vec<(u8, u8)> = Vec::new();
    let mut md_iter = md.iter().peekable();

    // TODO: clean this up, I think we can get rid of this outer while loop
    while let Some(next_char) = md_iter.next() {
        // This should be a number, it represents the number of matches
        let mut matches = vec![*next_char];

        // It could be more than a single number...
        while let Some(char_after) = md_iter.peek() {
            // Looks like the run number is more than a single digit
            if **char_after >= b'0' && **char_after <= b'9' {
                matches.push(**char_after);

                // Move to the next character
                md_iter.next();
            } else {
                break;
            }
        }

        // Add the number of matches to the expanded MD string
        for _ in 0..util::bytes_to_string(&matches).parse::<u32>().unwrap() {
            new_md.push((b'M', b'M'));
        }

        // The first set of characters we processed should have all been matches, now we
        // check for mismatches or deletions
        match md_iter.next() {
            // We're at the end of the MD tag, so we can break
            None => break,
            // This is a deletion, the deletion character should be followed by a string of bases
            Some(b'^') => {
                while let Some(del_char) = md_iter.peek() {
                    if **del_char >= b'A' && **del_char <= b'Z' {
                        new_md.push((b'D', **del_char));
                        md_iter.next();
                    } else {
                        break;
                    }
                }
            }
            // Either mismatches or we're back to a number which is signifies matches
            Some(char_after) => {
                // This is a mismatch
                if *char_after >= b'A' && *char_after <= b'Z' {
                    new_md.push((b'N', *char_after));

                    // continue plucking mismatches from the iterator
                    while let Some(mis_char) = md_iter.peek() {
                        if **mis_char >= b'A' && **mis_char <= b'Z' {
                            new_md.push((b'N', **mis_char));
                            md_iter.next();
                        } else {
                            break;
                        }
                    }
                }
            }
        }
    }

    new_md
    // Add any deletions or mismatches
}

/**
 * This uses the CIGAR string, MD tag, and read (query) sequence to reconstruct the
 * reference sequence and its alignment to the query sequence.
 *
 * TODO: this currently doesn't handle all CIGAR string characters.
 *
 * args
 *  cigar:    an expanded CIGAR string
 *  md:       an expanded MD tag
 *  sequence: the read (query) sequence
 *
 * returns
 *  a tuple containing
 *      0: the aligned and reconstructed reference sequence
 *      1: the aligned query sequence
 * both returned sequences are aligned to one another and should be equal length
 */
pub fn reconstruct_alignment(cigar: &[u8], md: &[(u8, u8)], sequence: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut reference = Vec::new();
    let mut query = Vec::new();
    let mut cigar_iter = cigar.iter().enumerate().peekable();
    let mut md_iter = md.iter().peekable();
    let mut seq_iter = sequence.iter().peekable();

    while cigar_iter.peek().is_some() && md_iter.peek().is_some() {
        let (cigar_ndx, cigar_char) = cigar_iter.peek().unwrap();
        let md_pair = **md_iter.peek().unwrap();

        match **cigar_char {
            // Deletion from reference
            b'D' => {
                reference.push(md_pair.1);
                query.push(b'-');

                md_iter.next();
                cigar_iter.next();
            }
            // Hard mask, nothing to do
            b'H' => {
                cigar_iter.next();
                continue;
            }
            // Insertion into reference
            b'I' => {
                reference.push(b'-');
                query.push(*seq_iter.next().unwrap());

                cigar_iter.next();
            }
            // Match
            b'M' => {
                if md_pair.0 == b'M' {
                    reference.push(**seq_iter.peek().unwrap());
                } else {
                    reference.push(md_pair.1);
                }

                query.push(*seq_iter.next().unwrap());

                cigar_iter.next();
                md_iter.next();
            }
            // Soft mask
            b'S' => {
                if md_pair.0 == b'M' {
                    reference.push(**seq_iter.peek().unwrap());
                } else {
                    reference.push(md_pair.1);
                }

                query.push(*seq_iter.next().unwrap());

                cigar_iter.next();
                md_iter.next();
            }
            x => panic!("{}", format!("Unhandled CIGAR Op: {}", x)),
        }
    }

    (reference, query)
}

/**
 * Kmerizes both the aligned reference and query sequences to produce a mapping of
 * every reference kmer to alternative kmers (determined via sequencing) and counts of those
 * alt kmers.
 *
 * args
 *  k:                kmer length
 *  alignment_length: length of the alignment
 *  reference:        reference sequence
 *  query:            query (read) sequence
 *
 * returns
 *  a mapping
 *      reference kmer -> alternate, aligned kmers and their counts
 */
pub fn encoded_kmerize_alignment(
    k: usize,
    alignment_length: usize,
    reference: &[u8],
    query: &[u8],
) -> EncodedKmerCountMap {
    let mut ndx: usize = 0;
    let mut kmer_map: EncodedKmerCountMap = HashMap::new();
    let mut invalid_ref_kmers: u32 = 0;
    let mut invalid_query_kmers: u32 = 0;

    while (ndx + k) < alignment_length {
        // If the reference kmer starts with a gap, we'll just skip and move on to the next
        if reference[ndx] == b'-' {
            ndx += 1;
            continue;
        }

        // Grab a kmer from the reference and query, remove any gaps along the way
        let ref_kmer = reference[ndx..(ndx + k)]
            .iter()
            .filter(|n| **n != b'-')
            .filter(|n| **n == b'A' || **n == b'C' || **n == b'G' || **n == b'T')
            .cloned()
            .collect::<Vec<u8>>();

        // For the query kmer we also remove gaps and any N bases
        let mut query_kmer = query[ndx..(ndx + k)]
            .iter()
            .filter(|n| **n != b'-' && **n != b'N')
            .cloned()
            .collect::<Vec<u8>>();

        // Skip reference kmers with gaps in them
        if ref_kmer.len() != k {
            ndx += 1;
            continue;
        }

        // Move to the next kmer if the query sequence is empty
        if query_kmer.is_empty() {
            ndx += 1;
            continue;
        }

        // Encode the kmers, if encoding fails due to invalid bases, e.g., currently this
        // doesn't support things IUPAC bases, then we just move onto the next kmer
        let ref_kmer_encoded = match encoding::three_bit_encode_kmer(&ref_kmer, k) {
            Ok(u) => u,
            Err(e) => {
                ndx += 1;
                continue;
            }
        };

        // Pad so the kmer is the same length as all others
        while query_kmer.len() < k {
            query_kmer.push(b'N');
        }

        let query_kmer_encoded = match encoding::three_bit_encode_kmer(&query_kmer, k) {
            Ok(u) => u,
            Err(e) => {
                ndx += 1;
                continue;
            }
        };

        // Add the ref kmer itself and its alt to the list of alts
        let mut kmer_counts = kmer_map.entry(ref_kmer_encoded).or_default();

        // Update counts for the reference kmer
        *kmer_counts.entry(ref_kmer_encoded).or_default() += 1;

        // then the query (sequenced) kmer
        *kmer_counts.entry(query_kmer_encoded).or_default() += 1;

        ndx += 1;
    }

    kmer_map
}

/**
 * Uses bincode to serialize a custom model and write it to the given output.
 */
//pub fn serialize_alignments_to_path(filepath: &Path, align: &AlignmentRecord) -> Result<(), String> {
pub fn serialize_alignments_to_path(
    filepath: &Path,
    align: Vec<AlignmentRecord>,
) -> Result<(), String> {
    // Serialize into a binary format
    let bytes = bincode::serialize(&align).unwrap();

    // Set up file for writing
    let mut file = match fs::OpenOptions::new()
        .write(true)
        .create(true)
        .open(filepath)
    {
        Ok(f) => f,
        Err(e) => return Err(format!("{}", e)),
    };

    // Write out the bytes to the file
    match file.write_all(&bytes) {
        Ok(_) => Ok(()),
        Err(e) => Err(format!("{}", e)),
    }
}

/**
 * Deserialize a custom model into the proper struct.
 */
pub fn deserialize_alignments_from_path(filepath: &Path) -> Result<Vec<AlignmentRecord>, String> {
    let bytes = match fs::read(filepath) {
        Ok(bs) => bs,
        Err(e) => return Err(format!("{}", e)),
    };

    // Attempt to deserialize
    let model: Vec<AlignmentRecord> = match bincode::deserialize(&bytes) {
        Ok(m) => m,
        Err(e) => return Err(format!("{}", *e)),
    };

    Ok(model)
}

pub fn serialize_to_hex_string<T: Serialize>(t: &T) -> Result<String, String> {
    let bytes = match bincode::serialize(&t) {
        Ok(b) => b,
        Err(e) => return Err(format!("{}", e)),
    };

    Ok(hex::encode(bytes))
}

pub fn deserialize_from_hex_string<T: DeserializeOwned>(s: &str) -> Result<T, String> {
    let bytes = match hex::decode(s) {
        Ok(b) => b,
        Err(e) => return Err(format!("{}", e)),
    };

    let t: T = match bincode::deserialize(&bytes) {
        Ok(t) => t,
        Err(e) => return Err(format!("{}", e)),
    };

    Ok(t)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_expand_cigar() {
        assert!(std::str::from_utf8(&expand_cigar("5M".as_bytes())).unwrap() == "MMMMM");
        assert!(std::str::from_utf8(&expand_cigar("2M1I".as_bytes())).unwrap() == "MMI");
        assert!(
            std::str::from_utf8(&expand_cigar("3H1M2D1I2M".as_bytes())).unwrap() == "HHHMDDIMM"
        );
    }

    #[test]
    fn test_expand_md_tag() {
        let truth0 = vec![
            (b'M', b'M'),
            (b'N', b'A'),
            (b'N', b'C'),
            (b'M', b'M'),
            (b'N', b'T'),
            (b'N', b'C'),
        ];
        let truth1 = vec![
            (b'M', b'M'),
            (b'M', b'M'),
            (b'N', b'G'),
            (b'N', b'A'),
            (b'M', b'M'),
            (b'D', b'A'),
            (b'D', b'T'),
            (b'M', b'M'),
        ];
        assert!(expand_md_tag("1A0C1T0C".as_bytes())
            .iter()
            .zip(truth0)
            .all(|(a, b)| a.0 == b.0 && a.1 == b.1));

        assert!(expand_md_tag("2G0A1^AT1".as_bytes())
            .iter()
            .zip(truth1)
            .all(|(a, b)| a.0 == b.0 && a.1 == b.1));
    }

    #[test]
    fn test_reconstruct_alignment() {
        // Additional tests I should put in

        //let cigar = expand_cigar("36M".as_bytes());
        //let md = expand_md_tag("1A0C0C0C1T0C0T27".as_bytes());
        //let query = "CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG".as_bytes();

        //let cigar = expand_cigar("6M1I29M".as_bytes());
        //let md = expand_md_tag("0C1C0C1C0T0C27".as_bytes());
        //let query = "GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT".as_bytes();

        //let cigar = expand_cigar("9M9D27M".as_bytes());
        //let md = expand_md_tag("2G0A5^ATGATGTCA27".as_bytes());
        //let query = "AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC".as_bytes();

        let cigar = expand_cigar("2M1I7M6D26M".as_bytes());
        let md = expand_md_tag("3C3T1^GCTCAG26".as_bytes());
        let query = "AGTGATGGGAGGATGTCTCGTCTGTGAGTTACAGCA".as_bytes();

        let query_truth = "AGTGATGGGA------GGATGTCTCGTCTGTGAGTTACAGCA".as_bytes();
        let ref_truth = "AG-GCTGGTAGCTCAGGGATGTCTCGTCTGTGAGTTACAGCA".as_bytes();

        let (aligned_ref, aligned_query) = reconstruct_alignment(&cigar, &md, query);

        assert!(aligned_ref.iter().zip(ref_truth).all(|(a, b)| a == b));
        assert!(aligned_query.iter().zip(query_truth).all(|(a, b)| a == b));
    }
}
