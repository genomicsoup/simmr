mod alignment;
mod cli;
mod log;
mod probability;

//use bincode::{config, Decode, Encode};
use itertools::Itertools;
use noodles::sam;
use noodles::sam::record::data::field::tag;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::{Arc, Mutex};
use tracing::{error, info, warn};

use shared::encoding;
use shared::util;

fn kmerize_alignments_from_disk(k: usize, alignment_path: &Path) -> alignment::EncodedKmerCountMap {
    let file = File::open(alignment_path).unwrap();

    //let mut kmer_map: alignment::EncodedKmerCountMap = HashMap::new();
    let mut kmer_map: Arc<Mutex<alignment::EncodedKmerCountMap>> =
        Arc::new(Mutex::new(HashMap::new()));

    io::BufReader::new(file)
        .lines()
        .par_bridge()
        .for_each(|line| {
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
            let kmers =
                alignment::encoded_kmerize_alignment(k, ref_seq.len(), &ref_seq, &query_seq);

            // We need to consolidate kmers across alignments, so merge all of the hashmaps together
            //kmers.into_iter().for_each(|(k, mut v)| {
            kmers.into_iter().for_each(|(k, v)| {
                let mut kmer_counts = kmer_map.lock().unwrap();
                let kmer_counts = kmer_counts.entry(k).or_default();

                for (subkey, c) in v.into_iter() {
                    *kmer_counts.entry(subkey).or_default() += c;
                }
            });
        });

    let final_map = kmer_map.lock().unwrap().clone();

    final_map
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
        kmers.into_iter().for_each(|(k, v)| {
            let kmer_counts = kmer_map.entry(k).or_default();

            for (subkey, c) in v.into_iter() {
                *kmer_counts.entry(subkey).or_default() += c;
            }
        });
    }

    kmer_map
}

/*
 * Generate and store all the distributions and relevant data necessary for read simulation.
 */
fn generate_simulation_distributions(args: &cli::CliArgs) {
    // Load the alignments from disk
    let mut sam_reader = File::open(&args.sam_file[0])
        .map(BufReader::new)
        .map(sam::Reader::new)
        .unwrap();

    let sam_header = sam_reader.read_header().unwrap();

    let temp_alignment_path = Path::new(&args.temp_directory).join("alignments.bin");

    info!("Parsing {}", args.sam_file[0]);

    let mut alignments = Vec::new();
    // Per-base quality scores
    let mut qualities: HashMap<u32, Vec<u8>> = HashMap::new();
    // Insert sizes for each paired alignment
    let mut insert_sizes: Vec<f64> = Vec::new();
    // Read lengths
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

        //if (i % 250_000) == 0 && i > 0 {
        //    info!("Processed {} records", i);
        //}

        // Record the base call positions and their quality scores
        for (i, score) in record.quality_scores().as_ref().iter().enumerate() {
            qualities
                .entry(i as u32)
                .or_default()
                .push(u8::from(*score));
        }

        // Get the read's sequence from the alignment record
        let seq = record.sequence().to_string().as_bytes().to_vec();
        let mapq = record.mapping_quality();

        // If a sequence isn't provided, skip. Probably an alignment w/ MAPQ == 0
        if seq.len() == 0 {
            missing_sequence += 1;
            continue;
        }

        // Store the read length
        read_lengths.push(seq.len() as f64);

        // Skip unmapped reads, these are only used for quality distributions
        if record.flags().is_unmapped() {
            unmapped_read += 1;
            continue;
        }

        // Also skip alignments with MAPQ == 0
        if mapq.is_some() && mapq.unwrap().get() == 0 {
            bad_quality_alignments += 1;
            continue;
        }

        // Template length of zero usually indicates that the mate is unmapped
        if !args.single_reads
            && record.template_length().abs() == 0
            && record.flags().is_mate_unmapped()
        {
            unmapped_mate_pairs += 1;
            continue;
        }

        // MD tag is required
        let md_tag = match record.data().get(&tag::MISMATCHED_POSITIONS) {
            None => {
                warn!(
                    "Read ({}) alignment is missing the MD tag",
                    util::bytes_to_string(record.read_name().unwrap().as_ref())
                );
                continue;
            }
            Some(t) => t.as_str().unwrap(),
        };

        // There are some (a few thousand during my experimentation) alignments with enormous
        // insert sizes, skip these b/c they distort the distribution
        if !args.single_reads && record.template_length().abs() > 5000 {
            continue;
        }

        // It's possible to have negative insert sizes if the first read is mapped
        // to the reverse strand.
        insert_sizes.push(record.template_length().abs() as f64);

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
            cigar,
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

    let kmer_map = if !args.in_memory {
        kmerize_alignments_from_disk(args.k, &temp_alignment_path)
    } else {
        kmerize_alignments(args.k, alignments)
    };

    info!(
        "Generating kmer probabilities for {} reference kmers",
        kmer_map.len()
    );

    let mut kmer_probs = probability::make_encoded_kmer_probabilities(kmer_map);

    // limit alternate kmers to the N most prevalent (highest prob) kmers
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

    info!("Generating quality score PDFs");

    let binned = probability::create_quality_bins(qualities, args.bin_size);

    info!("Generating read length and insert size PDFs");

    // Attempt to deterimn if these are long or short reads
    let is_long = util::mean(&read_lengths) > 400.0;

    // Need to sort the lengths prior to figuring out the pdf stuff
    read_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let read_length_bins = probability::create_read_length_bins(&read_lengths);
    let insert_size_bins = if insert_sizes.is_empty() || is_long {
        None
    } else {
        Some(probability::create_insert_size_bins(&insert_sizes))
    };

    info!("Model parameters:");
    info!("  read type: {}", if is_long { "long" } else { "short" });
    info!("  quality bin size: {}", args.bin_size);
    info!("  read length bin size: {}", read_length_bins.bin_width);
    info!("  read length bins: {}", read_length_bins.num_bins);
    info!("  k-mer bit encoding: {}", 3);
    info!("  k-mer size: {}", args.k);
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
            insert_size_bins,
            read_length_mean: util::mean(&read_lengths),
            read_length_std: util::std_deviation(&read_lengths),
            read_length_bins,
            is_long,
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

fn main() {
    let args = cli::parse_cli_args();

    // Setup stderr logging
    log::setup_logging();

    // Setup threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    generate_simulation_distributions(&args);

    //if args.generate_samples.is_some() {
    //    let model_params = match encoding::deserialize_model_from_path(path::Path::new(
    //        args.generate_samples.as_ref().unwrap(),
    //    )) {
    //        Ok(model) => model,
    //        Err(e) => {
    //            eprintln!("Error parsing custom error profile: {}", e);
    //            std::process::exit(1);
    //        }
    //    };
    //    let bins = &model_params.binned_quality_density[13191];

    //    let dist = WeightedAliasIndex::new(bins.binned_density.clone()).unwrap();
    //    let mut rng = StdRng::from_entropy();
    //    let score_choices = (0..=70).map(|x| x).collect_vec();
    //    let scores = (0..10_000)
    //        .map(|_| score_choices[dist.sample(&mut rng)])
    //        .collect_vec();

    //    let sp = Path::new("rando-scores.txt");
    //    let mut spfile = File::create(&sp).unwrap();

    //    scores.iter().for_each(|x| {
    //        writeln!(spfile, "{}", x).unwrap();
    //    });
    //    std::process::exit(0);
    //}
}
