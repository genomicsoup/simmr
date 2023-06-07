mod alignment;
mod cli;
mod log;
mod probability;

//use bincode::{config, Decode, Encode};
use bincode;
use itertools::Itertools;
use noodles::sam::record::cigar::op::kind::Kind;
//use noodles::sam::record::data::field::Tag;
use noodles::sam::record::data::field::tag;
use noodles::sam::{self, record::cigar};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::hash::Hash;
use std::io::Write;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use tracing::{error, info, warn};

use shared::encoding;
use shared::util;

fn convert_expanded_md_tag(md: Vec<(u8, u8)>) -> Vec<(char, char)> {
    md.iter().map(|(t, n)| (*t as char, *n as char)).collect()
}

fn kmerize_alignments_from_disk(k: usize, alignment_path: &Path) -> alignment::EncodedKmerCountMap {
    let file = File::open(alignment_path).unwrap();
    let mut kmer_map: alignment::EncodedKmerCountMap = HashMap::new();

    for (i, line) in io::BufReader::new(file).lines().enumerate() {
        // Deserialize the alignment record
        let rec: alignment::AlignmentRecord =
            alignment::deserialize_from_hex_string(&line.unwrap())
                .expect("Failed to deserialize alignment record");

        // Reconstruct the query/reference alignment using the expanded CIGAR and MD tag strings
        let (ref_seq, query_seq) = alignment::reconstruct_alignment(
            &alignment::expand_cigar(&rec.cigar),
            &alignment::expand_md_tag(&rec.md_tag),
            &rec.sequence,
        );

        // Kmerize the alignment and generate a map of each kmer -> sequenced alternate kmers
        let kmers = alignment::encoded_kmerize_alignment(k, ref_seq.len(), &ref_seq, &query_seq);

        // We need to consolidate kmers across alignments, so merge all of the hashmaps together
        kmers.into_iter().for_each(|(k, mut v)| {
            let kmer_counts = kmer_map.entry(k).or_default();

            for (subkey, c) in v.into_iter() {
                *kmer_counts.entry(subkey).or_default() += c;
            }
        });
    }

    kmer_map
}

fn kmerize_alignments(
    k: usize,
    alignments: Vec<alignment::AlignmentRecord>,
) -> alignment::EncodedKmerCountMap {
    let mut kmer_map: alignment::EncodedKmerCountMap = HashMap::new();

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
        let kmers = alignment::encoded_kmerize_alignment(k, ref_seq.len(), &ref_seq, &query_seq);

        // We need to consolidate kmers across alignments, so merge all of the hashmaps together
        kmers.into_iter().for_each(|(k, mut v)| {
            let kmer_counts = kmer_map.entry(k).or_default();

            for (subkey, c) in v.into_iter() {
                *kmer_counts.entry(subkey).or_default() += c;
            }
        });
    }

    kmer_map
}

fn main() {
    let args = cli::parse_cli_args();

    // Setup stderr logging
    log::setup_logging();

    let mut sam_reader = File::open(&args.sam_file[0])
        .map(BufReader::new)
        .map(sam::Reader::new)
        .unwrap();

    //let sam_header = sam_reader.read_header().unwrap().parse().unwrap();
    let sam_header = sam_reader.read_header().unwrap();

    let temp_alignment_path = Path::new(&args.temp_directory).join("alignments.bin");

    info!("Parsing {}", args.sam_file[0]);

    let mut alignments = Vec::new();
    // Save per-base quality scores
    let mut qualities: HashMap<u32, Vec<u8>> = HashMap::new();
    // Save insert sizes for each paired alignment
    let mut insert_sizes: Vec<f64> = Vec::new();
    // Save read lengths
    let mut read_lengths: Vec<f64> = Vec::new();
    // Inform users of alignments where MAPQ is 0
    let mut bad_quality_alignments = 0;
    // Inform users where the mate pair is unmapped
    let mut unmapped_mate_pairs = 0;
    // Inform users when the read itself is unmapped
    let mut unmapped_read = 0;
    // Inform users when the sequence is missing
    let mut missing_sequence = 0;

    for (i, res) in sam_reader.records(&sam_header).enumerate() {
        let record = res.unwrap();

        // Stop collecting alignments if necessary
        if args.max_alignments.is_some() && i > args.max_alignments.unwrap() {
            break;
        }

        if (i % 250_000) == 0 && i > 0 {
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
            unmapped_read += 1;
            continue;
        }

        // Get the read's sequence from the alignment record
        let seq = record.sequence().to_string().as_bytes().to_vec();
        let mapq = record.mapping_quality();

        // Also skip alignments with MAPQ == 0
        if mapq.is_some() && mapq.unwrap().get() == 0 {
            bad_quality_alignments += 1;
            continue;
        }

        // If a sequence isn't provided, skip. Probably an alignment w/ MAPQ == 0
        if seq.len() == 0 {
            missing_sequence += 1;
            continue;
        }

        // Template length of zero usually indicates that the mate is unmapped, alignments
        // where the mate is unmapped
        if record.template_length().abs() == 0 && record.flags().is_mate_unmapped() {
            unmapped_mate_pairs += 1;
            continue;
        }

        // MD tag is required
        //let md_tag = match record.data().get(Tag::MismatchedPositions) {
        let md_tag = match record.data().get(&tag::MISMATCHED_POSITIONS) {
            None => {
                warn!(
                    "Read ({}) alignment is missing the MD tag",
                    util::bytes_to_string(record.read_name().unwrap().as_ref())
                );
                continue;
            }
            //Some(t) => t.value().as_str().unwrap(),
            Some(t) => t.as_str().unwrap(),
        };

        if record.template_length().abs() == 0 && !record.flags().is_mate_unmapped() {
            println!("{:?}", record.template_length().abs());
            println!("{:?}", mapq);
            println!("qc fail {}", record.flags().is_qc_fail());
            println!("unmapped {}", record.flags().is_unmapped());
            println!("mate unmapped {}", record.flags().is_mate_unmapped());
            println!("mate align start {:?}", record.mate_alignment_start());
            println!("align start {:?}", record.alignment_start());
            println!("is_secondary {:?}", record.flags().is_secondary());
            println!("is_segment {:?}", record.flags().is_segmented());
            std::process::exit(1);
        }

        // There are some (a few thousand during my experimentation) alignments with enourmous
        // insert sizes, skip these b/c they distort the distribution
        if record.template_length().abs() > 5000 {
            continue;
        }

        // It's possible to have negative insert sizes if the first read is mapped
        // to the reverse strand.
        insert_sizes.push(record.template_length().abs() as f64);
        read_lengths.push(seq.len() as f64);

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

    info!("Using {} alignments", alignments.len());
    info!(
        "Skipped {} alignments with MAPQ == 0",
        bad_quality_alignments
    );
    info!(
        "Skipped {} alignments that were missing sequences",
        missing_sequence
    );
    info!(
        "Skipped {} alignments where the read was unmapped",
        unmapped_read
    );
    info!(
        "Skipped {} alignments where the mate was unmapped",
        unmapped_mate_pairs
    );

    // We're not doing this completely in memory, so store alignments to disk
    if !args.in_memory {
        info!(
            "Writing alignments to temporary location, {}",
            temp_alignment_path.display()
        );

        let mut f = File::create(&temp_alignment_path).unwrap();

        // Write encoded alignments to disk
        alignments
            .iter()
            .map(|a| alignment::serialize_to_hex_string(&a).unwrap())
            .for_each(|a| writeln!(f, "{}", a).unwrap());

        alignments.clear();
    }

    info!("Kmerizing alignments and encoding kmers");

    let mut kmer_map = if !args.in_memory {
        kmerize_alignments_from_disk(args.k, &temp_alignment_path)
    } else {
        kmerize_alignments(args.k, alignments)
    };

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

    // write insert sizes to disk
    let insert_size_path = Path::new("sizes.txt");
    let mut insert_size_file = File::create(&insert_size_path).unwrap();
    insert_sizes
        .iter()
        .for_each(|s| writeln!(insert_size_file, "{}", s).unwrap());

    info!("Model parameters:");
    info!("  bin size: {}", args.bin_size);
    info!("  bit encoding: {}", 3);
    info!("  kmer size: {}", args.k);
    info!("  insert size mean: {}", util::mean(&insert_sizes));
    info!("  insert size std: {}", util::std_deviation(&insert_sizes));
    info!("  read length mean: {}", util::mean(&read_lengths));
    info!("  read length std: {}", util::std_deviation(&read_lengths));

    let output_result = encoding::serialize_model_to_path(
        Path::new(&args.output),
        &encoding::ErrorModelParams {
            bin_size: args.bin_size,
            binned_quality_density: binned,
            bit_encoding: 3,
            kmer_size: args.k,
            probabilities: kmer_probs,
            insert_size_mean: util::mean(&insert_sizes),
            insert_size_std: util::std_deviation(&insert_sizes),
            read_length_mean: util::mean(&read_lengths),
            read_length_std: util::std_deviation(&read_lengths),
            is_long: util::mean(&read_lengths) > 400.0,
        },
    );

    // Clean up any temp files
    if temp_alignment_path.exists() {
        fs::remove_file(&temp_alignment_path).unwrap();
    }

    info!("Wrote sequence error model to {}", args.output);

    if output_result.is_err() {
        error!("Failed to ")
    }
}
