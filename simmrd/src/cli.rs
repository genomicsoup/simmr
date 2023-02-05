/**
 * file: cli.rs
 * desc: CLI parsing.
 */
use clap::Parser;

/**
 * STRUCTS
 */

#[derive(Debug, Parser)]
#[clap(version, about, long_about = None)]
pub struct CliArgs {
    #[clap(required(true), long, value_parser, help = "SAM file")]
    pub sam_file: Vec<String>,

    #[clap(required(true), long, value_parser, help = "Output file")]
    pub output: String,

    #[clap(
        long,
        value_parser,
        default_value_t = 5,
        help = "Quality score bin size used for kernel density estimation"
    )]
    pub bin_size: usize,

    #[clap(
        long,
        value_parser,
        help = "MAPQ threshold, alignments below the threshold will not be used"
    )]
    pub mapq_threshold: Option<u8>,

    #[clap(
        long,
        value_parser,
        help = "Use a maximum of N alignments for distribution modeling"
    )]
    pub max_alignments: Option<usize>,

    #[clap(
        long,
        value_parser,
        default_value_t = 20,
        help = "Use a maximum of N alternately sequenced kmers for distribution modeling"
    )]
    pub max_alt_kmers: u8,

    #[clap(
        long,
        // This doesn't work with usize
        //value_parser = clap::value_parser!(usize).range(3..10),
        value_parser = valid_kmer_size,
        default_value_t = 7,
        help = "Kmer length, must be between 3 and 10"
    )]
    pub k: usize,
}

/**
 * FUNCTIONS
 */

pub fn parse_cli_args() -> CliArgs {
    CliArgs::parse()
}

fn valid_kmer_size(s: &str) -> Result<usize, String> {
    let kmer_size: usize = s
        .parse()
        .map_err(|_| format!("`{}` isn't a valid integer", s))?;

    if kmer_size >= 3 && kmer_size <= 10 {
        Ok(kmer_size as usize)
    } else {
        Err(format!("Kmer size must be between {}-{}", 3, 10))
    }
}
