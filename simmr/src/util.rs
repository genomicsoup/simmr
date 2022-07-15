/**
 * file: util.rs
 * desc: Misc. utility functions.
 */
use uuid::Uuid;

/**
 * sequence-related functions
 */

/**
 * Generate the complement of the given nucleotide.
 * TODO: IUPAC support
 */
pub fn complement(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        x => x,
    }
}

/**
 * Given a sequence of nucleotides, generate the reverse complement.
 */
pub fn reverse_complement(nucs: &[u8]) -> Vec<u8> {
    nucs.iter()
        .rev()
        .map(|n| complement(*n))
        .collect::<Vec<u8>>()
}

/**
 * quality score functions
 */

/**
 * Encode a phred quality score into its ascii equivalent.
 */
pub fn encode_quality_score(s: u8) -> u8 {
    const PHRED_OFFSET: u8 = 33;

    s + PHRED_OFFSET
}

/**
 * Encode an array of phred quality score into their ascii representations.
 */
pub fn encode_quality_scores(scores: &Vec<u8>) -> Vec<u8> {
    scores.iter().map(|q| encode_quality_score(*q)).collect()
}

/**
 * Convert a phred quality score to an error probability.
 * https://gatk.broadinstitute.org/hc/en-us/articles/360035531872-Phred-scaled-quality-scores
 *
 * args
 *  score: a phred quality score
 *
 * returns
 *  an error probability (miscalled base)
 */
pub fn convert_phred_to_probability(score: u8) -> f32 {
    (10.0 as f32).powf(-((score as f32) / 10.0))
}

/**
 * Convert an error probability to a phred quality score.
 *
 * args
 *  prob: an error probability
 *
 * returns
 *  a phred quality score
 */
#[allow(dead_code)]
pub fn convert_probability_to_phred(prob: f32) -> u8 {
    (-10.0 * (prob).log10()) as u8
}

/**
 * Convert a phred quality score to the accuracy (1 - error) of the called base.
 *
 * args
 *  score: a phred quality score
 *
 * returns
 *  the accurarcy
 */
pub fn convert_phred_to_accuracy(score: u8) -> f32 {
    1.0 - convert_phred_to_probability(score)
}

/**
 * Convert an accuracy to a phred quality score.
 *
 * args
 *  acc: accuracy
 *
 * returns
 *  phred quality score
 */
pub fn convert_accuracy_to_phred(acc: f32) -> u8 {
    (-10.0 * (1.0 - acc).log10()).round() as u8
}

/**
 * misc functions
 */

/**
 * Generate a random UUID. The UUID generated is the first 64 most significant bits of a
 * UUID v4. An optional set of bytes can be provided to make use of a fixed random seed.
 *
 * returns
 *  a 64 bit unique ID
 */
pub fn generate_id(bytes: Option<&[u8; 16]>) -> u64 {
    match bytes {
        Some(bs) => Uuid::from_bytes(*bs).as_u64_pair().0,
        _ => Uuid::new_v4().as_u64_pair().0,
    }
}

#[cfg(test)]
#[path = "tests/util_tests.rs"]
mod util_tests;
