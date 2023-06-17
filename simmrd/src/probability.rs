/**
 * file: alignment.rs
 * desc: Functions for calculating kmer/alignment probabilities and distributions.
 */
use num_traits::Num;
use rayon::prelude::*;
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

/**
 * Calculate the IQR for a set of data. The data must be sorted prior to using this function.
 */
pub fn interquartile_range<T>(data: &[T]) -> f64
where
    T: Num + Copy + Into<f64>,
{
    let q1 = data[(data.len() as f64 * 0.25).floor() as usize];
    let q3 = data[(data.len() as f64 * 0.75).floor() as usize];

    q3.into() - q1.into()
}

pub fn freedman_diaconis_rule<T>(data: &[T]) -> usize
where
    T: Num + Copy + Into<f64>,
{
    let iqr = interquartile_range(data);
    let n = data.len() as f64;

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
pub fn create_quality_bins(quals: &HashMap<u32, Vec<u8>>, bin_size: usize) -> Vec<encoding::Bins> {
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

    // Generate the range of values in each bin, for quality scores we just make one bin per score
    let bin_ranges = (0..num_bins)
        .collect::<Vec<usize>>()
        .iter_mut()
        .enumerate()
        .map(|(i, _)| {
            //let start = (i * bin_size) as u32;
            //let end = (((i + 1) * bin_size) as u32) - 1;

            //(start, end)
            (i as u32, i as u32)
        })
        .collect::<Vec<(u32, u32)>>();
    //println!("{:?}", bin_ranges.len());
    //println!("bin ranges: {:?}", bin_ranges);

    quals.into_iter().for_each(|(ndx, scores)| {
        let float_scores = scores.iter().map(|v| *v as f64).collect::<Vec<f64>>();
        // Estimate bandwidth
        let bandwidth = calculate_bandwidth(&float_scores);

        // Estimate density for each quality score using KDE, assuming a normal dist., and assuming the
        // only quality scores are 0 - 70
        let density = (0..=MAX_PHRED_SCORE)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|s| gaussian(*s as f64, &float_scores, bandwidth))
            .collect::<Vec<f64>>();

        bins[*ndx as usize] = encoding::Bins {
            num_bins,
            bin_width: bin_size,
            binned_density: density,
            bin_ranges: bin_ranges.clone(),
        }
    });

    bins
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
pub fn create_read_length_bins(lengths: &[f64]) -> encoding::Bins {
    let bin_size = match freedman_diaconis_rule(&lengths) {
        bs if bs > 1 => bs,
        _ => 10,
    };
    let min_length = lengths.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_length = lengths.iter().fold(f64::MIN, |a, &b| a.max(b));
    let num_bins = ((max_length - min_length) / bin_size as f64).ceil() as usize;

    let mut bins = Vec::<encoding::Bins>::with_capacity(num_bins as usize);
    bins.resize_with(num_bins, Default::default);

    // Generate the range of values in each bin
    let bin_ranges = (0..num_bins)
        .collect::<Vec<usize>>()
        .iter_mut()
        .enumerate()
        .map(|(i, _)| {
            let start = min_length as u32 + (i * bin_size) as u32;
            let end = min_length as u32 + (((i + 1) * bin_size) as u32) - 1;

            (
                start,
                if end > max_length as u32 {
                    max_length as u32
                } else {
                    end
                },
            )
        })
        .collect::<Vec<(u32, u32)>>();
    //println!("min length: {}", min_length);
    //println!("max length: {}", max_length);
    //println!("bin range length: {:?}", bin_ranges.len());
    //println!("bin ranges: {:?}", bin_ranges);

    // Estimate bandwidth
    let bandwidth = calculate_bandwidth(&lengths);

    // Estimate density for each the midpoint of each read length bin using KDE, assuming
    // a normal dist., etc
    let densities = bin_ranges
        .par_iter()
        .map(|(s, e)| gaussian((s + e) as f64 / 2.0, &lengths, bandwidth))
        .collect::<Vec<f64>>();

    encoding::Bins {
        num_bins,
        bin_width: bin_size,
        binned_density: densities,
        bin_ranges: bin_ranges.clone(),
    }
}

/**
 * Use KDE to estimate a probability density function for insert sizes.
 * Generates a vector of densities for bins of insert sizes.
 *
 * args
 *  quals:    a mapping of position -> quality score list (all qualities observed at that position)
 *  bin_size: the size of each quality bin
 *
 * returns
 *  a vector of quality bins. Each index in the outer vector is the base pair position, and the
 *  index of each inner vector is the bin number.
 */
pub fn create_insert_size_bins(sizes: &[f64]) -> encoding::Bins {
    let bin_size = match freedman_diaconis_rule(&sizes) {
        bs if bs > 1 => bs,
        _ => 10,
    };
    let min_length = sizes.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max_length = sizes.iter().fold(f64::MIN, |a, &b| a.max(b));
    let num_bins = ((max_length - min_length) / bin_size as f64).ceil() as usize;

    let mut bins = Vec::<encoding::Bins>::with_capacity(num_bins as usize);
    bins.resize_with(num_bins, Default::default);

    // Generate the range of values in each bin
    let bin_ranges = (0..num_bins)
        .collect::<Vec<usize>>()
        .iter_mut()
        .enumerate()
        .map(|(i, _)| {
            let start = min_length as u32 + (i * bin_size) as u32;
            let end = min_length as u32 + (((i + 1) * bin_size) as u32) - 1;

            (start, end)
        })
        .collect::<Vec<(u32, u32)>>();

    // Estimate bandwidth
    let bandwidth = calculate_bandwidth(&sizes);

    // Estimate density for each the midpoint of each read length bin using KDE, assuming
    // a normal dist., etc
    let densities = bin_ranges
        .par_iter()
        .map(|(s, e)| gaussian((s + e) as f64 / 2.0, &sizes, bandwidth))
        .collect::<Vec<f64>>();

    encoding::Bins {
        num_bins,
        bin_width: bin_size,
        binned_density: densities,
        bin_ranges: bin_ranges.clone(),
    }
}
