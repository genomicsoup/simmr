/**
 * file: main.rs
 * desc: it's the main
 */
mod abundance_profiles;
mod cli;
mod error_profiles;
mod fastq;
mod files;
mod genome;
mod log;
mod simulate;
mod util;

use std::path::Path;
use tracing::{error, info, warn};

use crate::error_profiles::ErrorProfile;

fn run_main() {
    let args = cli::parse_cli_args();

    // Set up logging
    log::setup_logging();

    info!("Loading genomes");

    // Parse genomes from cmd args or the genome file
    let mut genomes: Vec<genome::Genome> = if args.genome_file.is_some() {
        // Read genomes + metadata from the genome file
        let records = match files::parse_genome_file(Path::new(&args.genome_file.as_ref().unwrap()))
        {
            Ok(rs) => rs,
            Err(e) => {
                error!("Failed to read genome file: {}", e);
                std::process::exit(1);
            }
        };

        // Ensure the genome filepaths actually exist
        records.iter().for_each(|rec| {
            if !Path::new(&rec.filepath).exists() {
                error!("Genome ({}) does not exist", rec.filepath);
                std::process::exit(1);
            }
        });

        // Parse genomes
        records
            .iter()
            .map(|rec| {
                // Try to parse
                let result = genome::Genome::from_fasta(&rec.filepath);

                // Fail on any errors
                if result.is_err() {
                    error!("Failed to parse {}: {}", rec.filepath, result.unwrap_err());
                    std::process::exit(1);
                }

                let mut genome = result.unwrap();

                // Update the UUID to match a provided one if necessary
                if rec.uuid.is_some() {
                    genome.uuid = genome::UUID::from(rec.uuid.as_ref().unwrap().clone());
                }

                genome
            })
            .collect()
    } else {
        // Read+parse genomes provided via the cmd line args
        args.genome
            .iter()
            .map(|g| {
                // Try to parse
                let result = genome::Genome::from_fasta(&g);

                // Fail on any errors
                if result.is_err() {
                    error!("Failed to parse {}: {}", g, result.unwrap_err());
                    std::process::exit(1);
                }

                result.unwrap()
            })
            .collect()
    };

    // Figure out which error and abundance profiles to use
    let eprofile = cli::determine_error_profile(&args);
    let aprofile = cli::determine_abundance_profile(&args);

    info!("Ensuring genomes meet minimum sequence length requirements for simulation");

    // Exclude any sequences which don't meet the required minimum size for the given error profile
    genomes.iter_mut().for_each(|g| {
        g.sequence = g
            .sequence
            .iter()
            .filter(|s| {
                if s.size <= eprofile.minimum_genome_size().into() {
                    warn!(
                        "({}) Sequence {} doesn't meet size requirements, size = {}, min size = {}",
                        g,
                        std::str::from_utf8(&s.id).unwrap(),
                        s.size,
                        eprofile.minimum_genome_size()
                    );
                    false
                } else {
                    true
                }
            })
            .cloned()
            .collect();
    });

    // Now exclude genomes that no longer have sequences b/c they didn't meet the min length reqs
    genomes = genomes
        .iter()
        .filter(|g| {
            if g.sequence.len() == 0 {
                warn!(
                    "Removing {} from simulation, it doesn't have usable sequences",
                    g
                );
                false
            } else {
                true
            }
        })
        .cloned()
        .collect();

    // Update genome lengths in case things were filtered out
    genomes
        .iter_mut()
        .for_each(|g| g.num_seqs = g.sequence.len());

    info!("Simulating reads");

    // Simulate reads
    let reads = if eprofile.is_long_read() {
        simulate::simulate_long_reads(args.num_reads, &genomes, &eprofile, &aprofile, args.seed)
    } else {
        simulate::simulate_pe_reads(args.num_reads, &genomes, &eprofile, &aprofile, args.seed)
    };

    info!("Writing simulated reads to {}", args.output);

    // If the output already exists, delete it
    if Path::new(&args.output).exists() {
        std::fs::remove_file(&args.output).ok();
    }

    // If the metadata already exists delete it
    if Path::new(&format!("{}.tsv", args.output)).exists() {
        std::fs::remove_file(&format!("{}.tsv", args.output)).ok();
    }

    // Save reads to the output filepath
    reads.iter().for_each(|(_, uuid, _, _, rs)| {
        match fastq::write_to_fastq(uuid, rs, &args.output, &args.read_header_format, true) {
            Ok(o) => o,
            Err(e) => error!("Failed to write reads to the output file: {}", e),
        };
    });

    info!(
        "Writing simulation metadata to {}",
        format!("{}.tsv", args.output)
    );

    let metadata = reads
        .iter()
        .map(|(genome_path, uuid, num_reads, abundance, _)| {
            (
                uuid.clone(),
                genome_path.as_os_str().to_str().unwrap().to_string(),
                *num_reads,
                *abundance,
            )
        })
        .collect::<Vec<(genome::UUID, String, usize, f64)>>();

    // TODO: refactor this garbage
    /*
    let metadata = reads
        .iter()
        .map(|(genome_path, uuid, num_reads, abundance, _)| {
            // This is ugly and unoptimized but works for now
            // TODO: redo this
            let seq_ids: Vec<u64> = genomes
                .iter()
                .filter(|g| g.uuid == *uuid)
                .map(|g| g.sequence.iter().map(|s| s.uuid).collect::<Vec<u64>>())
                .flatten()
                .collect();

            // Associate sequence IDs with genome IDs and metadata
            seq_ids
                .iter()
                .map(|sid| {
                    (
                        *sid,
                        uuid.clone(),
                        genome_path.as_os_str().to_str().unwrap().to_string(),
                        *num_reads,
                        *abundance,
                    )
                })
                .collect::<Vec<(u64, genome::UUID, String, usize, f64)>>()
        })
        .flatten()
        .collect::<Vec<(u64, genome::UUID, String, usize, f64)>>();
        */

    // Write metadata to disk
    let meta_result = files::write_metadata(&metadata, &format!("{}.tsv", args.output));

    if meta_result.is_err() {
        error!(
            "Failed to write metadata file: {}",
            meta_result.unwrap_err()
        );
    }

    // Format for metadata
}

fn main() {
    run_main();
}
