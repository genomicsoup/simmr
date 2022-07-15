/**
 * file: alignment.rs
 * desc: Functions for calculating kmer/alignment probabilities and distributions.
 */
use std::collections::HashMap;

use crate::alignment::EncodedKmerCountMap;

/**
 * Convert the counts of reference -> alternate kmers into probabilities.
 *
 * args
 *  kmer_map: reference -> alt kmer map and counts
 *
 * returns
 *  a vector of alternate kmer probabilities
 */
pub fn make_encoded_kmer_probabilities(
    kmer_map: EncodedKmerCountMap,
) -> Vec<(u32, Vec<(u32, f32)>)> {
    let mut kmer_probs = Vec::new();

    // Generate reference -> alternate kmer probs.
    for (i, (refs, alts)) in kmer_map.into_iter().enumerate() {
        // Vectorize the counts
        let counts: Vec<(u32, u32)> = alts.into_iter().collect();
        let total = counts.iter().map(|t| t.1).sum::<u32>() as f32;

        // Convert to probabilities
        kmer_probs.push((
            refs,
            counts
                .into_iter()
                .map(|(k, c)| (k, c as f32 / total as f32))
                .collect(),
        ))
    }

    kmer_probs
}

/**
 * For a typical sequencing run, there are GBs worth of quality scores. Typically this too much
 * to save and use quality score density estimation. Instead, we bin the quality scores and save
 * the number of elements in each bin, then use those bins/counts for KDE. This function bins
 * the quality scores and returns those bins/counts.
 *
 * args
 *  quals:    a mapping of position -> quality score list (all qualities observed at that position)
 *  bin_size: the size of each quality bin
 *
 * returns
 *  a vector of quality bins. Each index in the outer vector is the base pair position, and the
 *  index of each inner vector is the bin number.
 */
pub fn create_quality_bins(quals: HashMap<u32, Vec<u8>>, bin_size: usize) -> Vec<Vec<u32>> {
    // Assuming the max quality score possible is 70, bins are based around that. In reality, the
    // score could be higher (i've seen as high as 90). For those scores, we just lump them into
    // the final bin.
    const MAX_PHRED_SCORE: usize = 70;

    let num_bins = (MAX_PHRED_SCORE as f32 / bin_size as f32).ceil() as usize;

    // Vector indices correspond base pair position, elements are bin counts
    let mut qual_bins: Vec<Vec<u32>> =
        Vec::with_capacity((*quals.keys().max().unwrap() + 1) as usize);

    // apparently with_capacity doesn't actually change the vector length so do that here
    qual_bins.resize((*quals.keys().max().unwrap() + 1) as usize, vec![]);

    for (ndx, scores) in quals {
        // Generate the bins for this index. This
        let mut bins: Vec<u32> = vec![0; num_bins];

        // Since all our bins are equal width, this should put the score into the right bin
        for s in scores {
            let the_bin = (s as f32 / bin_size as f32).floor() as usize;

            // This only applies to scores higher than our MAX_PHRED_SCORE
            if the_bin >= num_bins {
                bins[num_bins - 1] += 1;
            } else {
                // We just need to count the number of elements in the bin
                bins[the_bin] += 1;
            }
        }

        qual_bins[ndx as usize] = bins;
    }

    qual_bins
}
