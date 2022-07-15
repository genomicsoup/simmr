/**
 * file: files.rs
 * desc: Functions related to file reading, writing, and parsing.
 */
use csv;
use serde::Deserialize;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use crate::genome;
//use std::fs;
//use std::io::Write;
//
//use crate::simulate::PERead;
//use crate::util;

#[derive(Debug, Deserialize)]
pub struct GenomeRecord {
    #[serde(alias = "path")]
    pub filepath: String,
    pub uuid: Option<String>,
    pub abundance: Option<f32>,
}

/**
 * Attempt to detect if the given genome file is the simple or complex variant.
 * This just reads a single line and checks if there's more than one column.
 */
fn is_simple_variant(filepath: &Path) -> bool {
    let reader = BufReader::new(File::open(filepath).expect("Failed to open genome file"));

    for line in reader.lines() {
        // Can we split this using the tab delimiter? If so, it's not the simple variant type
        if line.unwrap().split("\t").collect::<String>().len() > 1 {
            return false;
        }

        break;
    }

    true
}
/**
 * Parse the genomee file. The format is TSV and can be one of two variants: simple or complex.
 * A simple genome file is just a list of filepaths to genome FASTAs.
 * The complex format has a header and includes additional genome metadata. The header can have
 * the following fields, everything except the path is optional:
 *      path:      path to the genome FASTA
 *      uuid:      unique ID for this genome
 *      abundance: generates reads at this abundance level for this genome
 *
 */
pub fn parse_genome_file(filepath: &Path) -> Result<Vec<GenomeRecord>, String> {
    if !filepath.exists() {
        return Err(format!("Genome file does not exist"));
    }

    let reader_result = if is_simple_variant(filepath) {
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(filepath)
    } else {
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(filepath)
    };

    let mut reader = match reader_result {
        Ok(r) => r,
        Err(e) => return Err(format!("{}", e)),
    };
    let mut recs = Vec::new();

    // Deserialize into GenomeRecords
    for rec_result in reader.deserialize() {
        let rec: GenomeRecord = rec_result.unwrap();

        recs.push(rec);
    }

    Ok(recs)
}

/**
 * Write simulation metadata to the output file.
 *
 * args
 *  metadata: a vector of tuples containing the simulation metadata
 *      0: sequence ID the read originated from
 *      1: genome ID the read originated from
 *      2: filepath to the original genome
 *      3: number of reads generated for this genome
 *      4: abundance of this genome in the simulated read dataset
 */

pub fn write_metadata(
    //metadata: &Vec<(u64, genome::UUID, String, usize, f64)>,
    metadata: &Vec<(genome::UUID, String, usize, f64)>,
    output: &str,
) -> Result<(), String> {
    // Delete the metadata file if it already exists
    if Path::new(output).exists() {
        std::fs::remove_file(output).ok();
    }

    // Set up file for writing
    let mut file = match fs::OpenOptions::new().write(true).create(true).open(output) {
        Ok(f) => f,
        Err(e) => return Err(format!("{}", e)),
    };

    // Generate the header. Convert results into Options b/c I don't really care if the writes fail
    file.write_all(
        format!(
            "{}\t{}\t{}\t{}",
            "genome_id", "filepath", "num_reads", "abundance"
        )
        .as_bytes(),
    )
    .ok();
    file.write_all("\n".as_bytes()).ok();

    metadata.iter().for_each(|(gid, gpath, nreads, abundance)| {
        file.write_all(format!("{}\t{}\t{}\t{}", gid, gpath, nreads, abundance).as_bytes())
            .ok();
        file.write_all("\n".as_bytes()).ok();
    });

    Ok(())
}

/*
pub fn write_to_fastq(reads: &Vec<PERead>, output: &str, append: bool) -> Result<(), String> {
    // Set up file for writing
    let mut file = match fs::OpenOptions::new()
        .write(true)
        .append(append)
        .create(true)
        .open(output)
    {
        Ok(f) => f,
        Err(e) => return Err(format!("{}", e)),
    };

    for read in reads.iter() {
        let fwd_header = format!("@{}/1", std::str::from_utf8(&read.id).unwrap());
        let rev_header = format!("@{}/2", std::str::from_utf8(&read.id).unwrap());

        // write the forward read
        file.write_all(fwd_header.as_bytes());
        file.write_all("\n".as_bytes());
        file.write_all(&read.forward.sequence);
        file.write_all("\n".as_bytes());
        file.write_all("+".as_bytes());
        file.write_all("\n".as_bytes());
        file.write_all(&util::encode_quality_scores(&read.forward.quality));
        file.write_all("\n".as_bytes());

        // write the reverse read
        file.write_all(rev_header.as_bytes());
        file.write_all("\n".as_bytes());
        file.write_all(&read.reverse.sequence);
        file.write_all("\n".as_bytes());
        file.write_all("+".as_bytes());
        file.write_all("\n".as_bytes());
        file.write_all(&util::encode_quality_scores(&read.reverse.quality));
        file.write_all("\n".as_bytes());
    }

    Ok(())
} */
