/**
 * file: alignment.rs
 * desc: Functions for calculating kmer/alignment probabilities and distributions.
 */
use std::collections::HashMap;
use std::f64::consts::PI;

use crate::alignment::EncodedKmerCountMap;
use crate::encoding;
use crate::util;

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

pub fn interquartile_range(data: &[u8]) -> f32 {
    //let mut sorted = data.to_vec();
    //sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    //let q1 = sorted[(data.len() as f32 * 0.25).floor() as usize];
    //let q3 = sorted[(data.len() as f32 * 0.75).floor() as usize];
    let q1 = data[(data.len() as f32 * 0.25).floor() as usize];
    let q3 = data[(data.len() as f32 * 0.75).floor() as usize];
    println!("q1: {}, q3: {}", q1, q3);

    q3 as f32 - q1 as f32
}

pub fn freedman_diaconis_rule(data: &[u8]) -> usize {
    let iqr = interquartile_range(data);
    println!("iqr: {}", iqr);
    let n = data.len() as f32;

    (2.0 * (iqr / n.powf(1.0 / 3.0))) as usize
}

pub fn scotts_rule(data: &[u8]) -> usize {
    let std_dev = util::std_deviation(&data.iter().map(|v| *v as f64).collect::<Vec<f64>>());
    println!("std dev: {}", std_dev);
    let n = data.len() as f64;

    ((3.49 * std_dev) / n.powf(1.0 / 3.0)) as usize
}

// TODO: put this elsewhere
// Use the normal kernel to generate a value from the density function using the given value (x).
// Implements PDF formula here: https://en.wikipedia.org/wiki/Gaussian_function
fn gaussian(x: f64, xs: &[f64], bandwidth: f64) -> f64 {
    let sum: f64 = xs
        .iter()
        .map(|f| (x - f) / bandwidth)
        .map(|f| (-0.5 * f.powi(2)).exp())
        .sum();

    sum / ((2.0 * PI).sqrt() * xs.len() as f64 * bandwidth)
}

fn calculate_bandwidth(values: &[f64]) -> f64 {
    // Estimate the bandwidth parameter using Silverman's rule of thumb, not the best but it's fast
    let bandwidth = 1.06
        //* util::std_deviation(&counts.iter().map(|v| *v as f64).collect::<Vec<f64>>())
        * util::std_deviation(&values)
        * (values.len() as f64).powf(-1.0 / 5.0);

    bandwidth
}

/*
pub struct Bins {
    pub num_bins: usize,
    pub bin_width: u8,
    pub binned_density: Vec<f64>,
    pub bin_ranges: Vec<(u8, u8)>,
}
    */

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
pub fn create_quality_bins(quals: HashMap<u32, Vec<u8>>, bin_size: usize) -> Vec<encoding::Bins> {
    // Assuming the max quality score possible is 70, bins are based around that. In reality, the
    // score could be higher (i've seen as high as 90). For those scores, we just lump them into
    // the final bin.
    const MAX_PHRED_SCORE: usize = 70;

    //let num_bins = (MAX_PHRED_SCORE as f32 / bin_size as f32).ceil() as usize;
    let num_bins = MAX_PHRED_SCORE;

    // Vector indices correspond base pair position, elements are bin counts
    let mut qual_bins: Vec<Vec<u32>> =
        Vec::with_capacity((*quals.keys().max().unwrap() + 1) as usize);

    // apparently with_capacity doesn't actually change the vector length so do that here
    qual_bins.resize((*quals.keys().max().unwrap() + 1) as usize, vec![]);

    let mut bins =
        Vec::<encoding::Bins>::with_capacity((*quals.keys().max().unwrap() + 1) as usize);
    bins.resize_with(
        (*quals.keys().max().unwrap() + 1) as usize,
        Default::default,
    );
    //bins.resize_with((*quals.keys().max().unwrap() + 1) as usize, vec![]);
    //let mut bins = Vec::with_capacity((*quals.keys().max().unwrap() + 1) as usize);
    // Generate the range of values in each bin
    let bin_ranges = Vec::<(u8, u8)>::with_capacity(num_bins)
        .iter_mut()
        .enumerate()
        .map(|(i, _)| {
            let start = (i * bin_size) as u8;
            let end = (((i + 1) * bin_size) as u8) - 1;

            (start, end)
        })
        .collect::<Vec<(u8, u8)>>();

    for (ndx, scores) in quals {
        let float_scores = scores.iter().map(|v| *v as f64).collect::<Vec<f64>>();
        // Estimate bandwidth
        let bandwidth = calculate_bandwidth(&float_scores);

        // Estimate density for each quality score using KDE, assuming a normal dist., and assuming the
        // only quality scores are 0 - 70
        let pdf = (0..=MAX_PHRED_SCORE)
            .map(|s| gaussian(s as f64, &float_scores, bandwidth))
            .collect::<Vec<f64>>();
        println!("ndx: {}, pdf: {:?}", ndx, pdf);
        //let pdf = scores
        //    .iter()
        //    .map(|s| gaussian(*s as f64, &scores.iter().map(|v| *v as f64).collect::<Vec<f64>>(), bandwidth))
        //    .collect::<Vec<f64>>();
        bins[ndx as usize] = encoding::Bins {
            num_bins,
            bin_width: bin_size as u8,
            binned_density: pdf,
            bin_ranges: bin_ranges.clone(),
        }
    }

    /*
    // Generate the range of values in each bin
    let bin_ranges = Vec::<(u8, u8)>::with_capacity(num_bins)
        .iter_mut()
        .enumerate()
        .map(|(i, _)| {
            let start = (i * bin_size) as u8;
            let end = (((i + 1) * bin_size) as u8) - 1;

            (start, end)
        })
        .collect::<Vec<(u8, u8)>>();

    let mut bins = Vec::with_capacity((*quals.keys().max().unwrap() + 1) as usize);

    for (ndx, scores) in quals {
        // Generate the bins for this index. This
        let mut bin_densities: Vec<u32> = vec![0; num_bins];

        for s in &scores {
            let the_bin = (*s as f32 / bin_size as f32).floor() as usize;

            // This only applies to scores higher than our MAX_PHRED_SCORE
            if the_bin >= num_bins {
                bin_densities[num_bins - 1] += 1;
            } else {
                // We just need to count the number of elements in the bin
                bin_densities[the_bin] += 1;
            }
        }

        println!("ndx: {}", ndx);

        // Create the bin for this position
        bins[ndx as usize] = encoding::Bins {
            num_bins,
            bin_width: bin_size as u8,
            // Convert into frequencies, then densities
            binned_density: bin_densities
                .into_iter()
                .map(|v| ((v as f64) / scores.len() as f64) / bin_size as f64)
                .collect(),
            bin_ranges: bin_ranges.clone(),
        }
    }
    */

    bins
    /*
    for (ndx, scores) in quals {
        // Generate the bins for this index. This
        let mut bins: Vec<u32> = vec![0; num_bins];

        let mut observations = scores.clone();

        observations.sort();

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
    */
}
