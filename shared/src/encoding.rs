/**
 * file: encoding.rs
 * desc: Functions related to DNA encoding and data serialization.
 */
use bincode;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::Write;
use std::path::Path;

/**
 * MACROS
 */

macro_rules! create_mask {
    ($k:expr) => {{
        0xffffffff & !(3 << $k)
    }};
}

/**
 * Generates a bitmask with cleared bits starting at bit k and length (number of cleared bits) n.
 * Returns that mask with cleared bits.
 */
macro_rules! create_mask2 {
    ($k:expr, $n:expr) => {{
        let mut mask = 0xffffffff;

        for i in 0..$n {
            mask &= !(1 << ($k + i));
        }

        mask
    }};
}

macro_rules! set_bits {
    ($num:expr, $region:expr, $value:expr) => {{
        ($num & create_mask!($region)) | ($value << $region)
    }};
}

/**
 * For the given integer (num), set the series of bits (n) starting at the kth bit to the
 * given value. Return the modified integer.
 */
macro_rules! set_bits2 {
    ($num:expr, $k:expr, $n:expr, $value:expr) => {{
        ($num & create_mask2!($k, $n)) | ($value << $k)
    }};
}

macro_rules! extract_bits {
    ($num:expr, $region:expr) => {{
        ((1 << 2) - 1) & ($num >> $region)
    }};
}

/**
 * For the given integer (num), extract n bits starting at the kth bit.
 * Return the integer representing extracted bits.
 */
macro_rules! extract_bits2 {
    ($num:expr, $n:expr, $k:expr) => {{
        ((1 << $n) - 1) & ($num >> $k)
    }};
}

/**
 * STRUCTS
 */

/**
 * Models the binning of quality scores for a single position in a read.
 *
 * fields
 *  num_bins:       the number of bins
 *  bin_width:      width/size of each bin
 *  binned_density: density of quality scores in each bin
 *  bin_ranges:     range of quality scores in each bin
 */
#[derive(Clone, Default, Deserialize, Serialize, PartialEq, Debug)]
pub struct Bins {
    pub num_bins: usize,
    pub bin_width: usize,
    pub binned_density: Vec<f64>,
    pub bin_ranges: Vec<(u32, u32)>,
}

/**
 * A struct with custom error model parameters and distributions. This is serialized using
 * bincode into an output that can be used by simmr to model quality scores and
 * sequencing errors.
 *
 * fields
 *  bin_size:      the bin size used to generate KDE functions quality score simulation
 *  qualities:     quality score bins and counts
 *  bit_encoding:  bit encoding size (either 2 or 3 but almost certainly 3) for encoded kmers
 *  kmer_size:     kmer length
 *  probabilities: a map of reference kmers to alternately sequenced kmers and their probabilities
 */
#[derive(Clone, Deserialize, Serialize, PartialEq, Debug)]
pub struct ErrorModelParams {
    pub bin_size: usize,
    //pub qualities: Vec<Vec<u32>>,
    pub binned_quality_density: Vec<Bins>,
    pub bit_encoding: u8,
    pub kmer_size: usize,
    pub probabilities: Vec<(u32, Vec<(u32, f32)>)>,
    pub insert_size_mean: f64,
    pub insert_size_std: f64,
    pub insert_size_bins: Option<Bins>,
    pub read_length_mean: f64,
    pub read_length_std: f64,
    pub read_length_bins: Bins,
    pub is_long: bool,
}

impl std::fmt::Display for ErrorModelParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();

        s.push_str("<Error Model Params>\n");
        s.push_str(&format!("bin_size: {}\n", self.bin_size));
        s.push_str(&format!("bit_encoding: {}\n", self.bit_encoding));
        s.push_str(&format!("kmer_size: {}\n", self.kmer_size));
        s.push_str(&format!("insert_size_mean: {}\n", self.insert_size_mean));
        s.push_str(&format!("insert_size_std: {}\n", self.insert_size_std));
        s.push_str(&format!("read_length_mean: {}\n", self.read_length_mean));
        s.push_str(&format!("read_length_std: {}\n", self.read_length_std));
        s.push_str(&format!("is_long: {}\n", self.is_long));

        write!(f, "{}", s)
    }
}

/**
 * FUNCTIONS
 */

/**
 * Encode the given kmer into an unsigned 32bit integer. This function assumes that
 * the kmer size (k) is small enough to fit into a u32. This function supports both
 * 2bit encoding of A, C, G, T and 3bit encoding of A, C, G, T, N. This encoding size
 * can be changed using the esize argument. Any other nucleotides will throw an error.
 * 2bit encoding: a = 0 (00); c = 1 (01); g = 2 (10); t = 3 (11)
 * 3bit encoding: a = 0 (000); c = 1 (001); g = 2 (010); t = 3 (011); n = (100)
 */
fn encode_kmer(kmer: &[u8], k: usize, esize: usize) -> Result<u32, String> {
    let mut encoded: u32 = 0;

    for i in 0..k {
        match kmer[i] {
            // A is encoded as 00, so we can just move on
            b'A' => continue,
            b'C' => encoded = set_bits2!(encoded, i * esize, esize, 1),
            b'G' => encoded = set_bits2!(encoded, i * esize, esize, 2),
            b'T' => encoded = set_bits2!(encoded, i * esize, esize, 3),
            // Requires 3-bit encoding
            b'N' => {
                encoded = {
                    if esize < 3 {
                        return Err(format!("Cannot encode the base 'N' using {} bits", esize));
                    } else {
                        set_bits2!(encoded, i * esize, esize, 4)
                    }
                }
            }
            unknown => {
                return Err(format!(
                    "Cannot encode the following base: {} ({})",
                    unknown, unknown as char
                ));
            }
        }
    }

    Ok(encoded)
}

/**
 * Decode a u32 encoded kmer. Like the encode_kmer function, this supports 2bit and 3bit DNA
 * encoding. This function doesn't do any kind of encoding checks (e.g., you provide a 3bit
 * encoded kmer but set esize=2). If skip_n is true
 */
fn decode_kmer(kmer: u32, k: usize, esize: usize, skip_n: bool) -> Result<Vec<u8>, String> {
    let mut decoded = Vec::new();

    for i in 0..k {
        //match extract_bits!(kmer, i * esize) {
        match extract_bits2!(kmer, esize, i * esize) {
            0 => decoded.push(b'A'),
            1 => decoded.push(b'C'),
            2 => decoded.push(b'G'),
            3 => decoded.push(b'T'),
            // Requires 3-bit encoding, if using 2-bit encoding we'll never get to these cases
            4 => {
                if skip_n {
                    continue;
                } else {
                    decoded.push(b'N');
                }
            }
            // We should never get here
            n => return Err(format!("Decoded an invalid value: {}", n)),
        }
    }

    Ok(decoded)
}

/**
 * 2bit encode the given kmer of size k.
 */
pub fn two_bit_encode_kmer(kmer: &[u8], k: usize) -> Result<u32, String> {
    encode_kmer(kmer, k, 2)
}

/**
 * 3bit encode the given kmer of size k.
 */
pub fn three_bit_encode_kmer(kmer: &[u8], k: usize) -> Result<u32, String> {
    encode_kmer(kmer, k, 3)
}

/**
 * 2bit decode the given kmer of size k.
 */
pub fn two_bit_decode_kmer(kmer: u32, k: usize) -> Result<Vec<u8>, String> {
    decode_kmer(kmer, k, 2, false)
}

/**
 * 3bit decode the given kmer of size k. If skip_n is true, will remove any N characters
 * frome the decoded kmer.
 */
pub fn three_bit_decode_kmer(kmer: u32, k: usize, skip_n: bool) -> Result<Vec<u8>, String> {
    decode_kmer(kmer, k, 3, skip_n)
}

/**
 * Uses bincode to serialize a custom model and write it to the given output.
 */
pub fn serialize_model_to_path(filepath: &Path, model: &ErrorModelParams) -> Result<(), String> {
    // Serialize into a binary format
    let bytes = bincode::serialize(&model).unwrap();

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
pub fn deserialize_model_from_path(filepath: &Path) -> Result<ErrorModelParams, String> {
    let bytes = match fs::read(filepath) {
        Ok(bs) => bs,
        Err(e) => return Err(format!("{}", e)),
    };

    // Attempt to deserialize
    let model: ErrorModelParams = match bincode::deserialize(&bytes) {
        Ok(m) => m,
        Err(e) => return Err(format!("{}", *e)),
    };

    Ok(model)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_2bit_encoding() {
        assert!(two_bit_encode_kmer(&"ACGT".as_bytes(), 4).unwrap() == 0x000000E4);
        assert!(two_bit_encode_kmer(&"AAAAA".as_bytes(), 5).unwrap() == 0x00000000);
        assert!(two_bit_encode_kmer(&"TTATC".as_bytes(), 5).unwrap() == 0x000001CF);
        assert!(two_bit_encode_kmer(&"GCGCATCT".as_bytes(), 8).unwrap() == 0x0000DC66);
    }

    #[test]
    fn test_2bit_decoding() {
        assert!(
            std::str::from_utf8(&two_bit_decode_kmer(0x000000E4, 4).unwrap()).unwrap()
                == "ACGT".to_string()
        );
        assert!(
            std::str::from_utf8(&two_bit_decode_kmer(0x00000000, 5).unwrap()).unwrap()
                == "AAAAA".to_string()
        );
        assert!(
            std::str::from_utf8(&two_bit_decode_kmer(0x000001CF, 5).unwrap()).unwrap()
                == "TTATC".to_string()
        );
        assert!(
            std::str::from_utf8(&two_bit_decode_kmer(0x0000DC66, 8).unwrap()).unwrap()
                == "GCGCATCT".to_string()
        );
    }
}
