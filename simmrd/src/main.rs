mod alignment;
mod cli;
mod log;
mod probability;

//use bincode::{config, Decode, Encode};
use bincode;
use itertools::Itertools;
use noodles::sam::record::cigar::op::kind::Kind;
use noodles::sam::record::data::field::Tag;
use noodles::sam::{self, record::cigar};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::hash::Hash;
use std::io::Write;
use std::path::Path;
use std::{fs::File, io::BufReader};
use tracing::{error, info, warn};

use shared::encoding;
use shared::util;

fn convert_expanded_md_tag(md: Vec<(u8, u8)>) -> Vec<(char, char)> {
    md.iter().map(|(t, n)| (*t as char, *n as char)).collect()
}
fn main() {
    let args = cli::parse_cli_args();

    // Setup stderr logging
    log::setup_logging();

    // "/Users/tim/code/simmr/align-tests/test-missing-seq.sam"
    // "/Volumes/minidata/SRR13059359.sam"
    // "/Users/tim/code/simmr/align-tests/test.sam"
    let mut sam_reader = File::open(&args.sam_file[0])
        .map(BufReader::new)
        .map(sam::Reader::new)
        .unwrap();

    let sam_header = sam_reader.read_header().unwrap().parse().unwrap();

    info!("Parsing {}", args.sam_file[0]);

    let mut alignments = Vec::new();
    // Save per-base quality scores
    let mut qualities: HashMap<u32, Vec<u8>> = HashMap::new();

    for (i, res) in sam_reader.records(&sam_header).enumerate() {
        let record = res.unwrap();

        // Stop collecting alignments if necessary
        if args.max_alignments.is_some() && i > args.max_alignments.unwrap() {
            break;
        }

        if (i % 10_000) == 0 && i > 0 {
            info!("Processed {} records", i);
        }

        // Record the base call positions and their quality scores
        for (i, score) in record.quality_scores().as_ref().iter().enumerate() {
            qualities
                .entry(i as u32)
                .or_default()
                .push(u8::from(*score));
        }

        // Skip unmapped reads, these are only used for quality distributions
        if record.flags().is_unmapped() {
            continue;
        }

        // MD tag is required
        let md_tag = match record.data().get(Tag::MismatchedPositions) {
            None => {
                warn!(
                    "Read ({}) alignment is missing the MD tag",
                    util::bytes_to_string(record.read_name().unwrap().as_ref())
                );
                continue;
            }
            Some(t) => t.value().as_str().unwrap(),
        };

        // Get the read's sequence from the alignment record
        let seq = record.sequence().to_string().as_bytes().to_vec();

        // If a sequence isn't provided, skip. Probably an alignment w/ MAPQ == 0
        if seq.len() == 0 {
            continue;
        }

        // Regenerate the raw CIGAR string from the alignment record
        let cigar: Vec<u8> = record
            .cigar()
            .to_vec()
            .iter()
            .map(|op| {
                format!("{}{}", op.len(), char::from(op.kind()))
                    .as_bytes()
                    .to_vec()
            })
            .flatten()
            .collect();

        alignments.push(alignment::AlignmentRecord {
            cigar: cigar,
            sequence: if record.flags().is_reverse_complemented() {
                util::reverse_complement(&seq)
            } else {
                seq
            },
            md_tag: md_tag.as_bytes().to_vec(),
        });
    }

    info!("Kmerizing alignments and encoding kmers");

    let mut kmers = Vec::new();

    for (i, rec) in alignments.iter().enumerate() {
        if (i % 20_000) == 0 && i > 0 {
            info!("Kmerized {} alignments", i);
        }

        // Reconstruct the query/reference alignment using the expanded CIGAR and MD tag strings
        let (ref_seq, query_seq) = alignment::reconstruct_alignment(
            &alignment::expand_cigar(&rec.cigar),
            &alignment::expand_md_tag(&rec.md_tag),
            &rec.sequence,
        );

        // Kmerize the alignment and generate a map of each kmer -> sequenced alternate kmers
        kmers.push(alignment::encoded_kmerize_alignment(
            args.k,
            ref_seq.len(),
            &ref_seq,
            &query_seq,
        ));
    }

    info!("Consolidating reference and alternate kmers");

    let mut kmer_map: alignment::EncodedKmerCountMap = HashMap::new();

    // We need to consolidate kmers across alignments, so merge all of the hashmaps together
    for km in kmers {
        km.into_iter().for_each(|(k, mut v)| {
            let kmer_counts = kmer_map.entry(k).or_default();

            for (subkey, c) in v.into_iter() {
                *kmer_counts.entry(subkey).or_default() += c;
            }
        })
    }

    info!(
        "Generating kmer probabilities for {} reference kmers",
        kmer_map.len()
    );

    let mut kmer_probs = probability::make_encoded_kmer_probabilities(kmer_map);

    // limit alternate kmers to the N prevalent (highest prob) kmers
    kmer_probs = kmer_probs
        .iter()
        .map(|(k, alts)| {
            (
                *k,
                alts.into_iter()
                    .sorted_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .map(|t| *t)
                    .rev()
                    .take(args.max_alt_kmers as usize)
                    .collect::<Vec<(u32, f32)>>(),
            )
        })
        .collect::<Vec<(u32, Vec<(u32, f32)>)>>();

    info!("Binning quality scores");

    let binned = probability::create_quality_bins(qualities, args.bin_size);

    let output_result = encoding::serialize_model_to_path(
        Path::new(&args.output),
        &encoding::ErrorModelParams {
            bin_size: args.bin_size,
            qualities: binned,
            bit_encoding: 3,
            kmer_size: args.k,
            probabilities: kmer_probs,
        },
    );

    if output_result.is_err() {
        error!("Failed to ")
    }
}
