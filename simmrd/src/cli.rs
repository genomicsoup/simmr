/**
 * file: cli.rs
 * desc: CLI parsing.
 */
use clap::{Parser, Subcommand};

/**
 * STRUCTS
 */

//#[derive(Debug, Parser, Subcommand)]
#[derive(Debug, Parser)]
pub struct SimulateCommand {
    #[clap(long, value_parser, help = "Distribution blob")]
    pub distribution: String,
    #[clap(
        long,
        value_parser,
        help = "Simulate insert sizes and save to the given file"
    )]
    pub insert_size: Option<String>,
}

#[derive(Debug, Parser)]
//#[clap(version, about, long_about = None)]
//pub struct CliArgs {
//#[derive(Debug, Parser, Subcommand)]
pub struct GenerateCommand {
    #[clap(long, value_parser, help = "SAM file")]
    pub sam_file: Vec<String>,

    #[clap(long, value_parser, help = "Output file")]
    pub output: String,

    #[clap(long, value_parser, help = "View distribution blob information")]
    pub view: Option<String>,

    #[clap(
        long,
        value_parser,
        help = "Generate random samples from the given distribution [unimplemented]"
    )]
    pub generate_samples: Option<String>,

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

    #[clap(
        long,
        value_parser,
        default_value_t = false,
        help = "Alignment contains single ended or long reads"
    )]
    pub single_reads: bool,

    #[clap(
        long,
        value_parser,
        default_value = "/tmp",
        help = "Temporary directory for intermediate files"
    )]
    pub temp_directory: String,

    #[clap(
        long,
        value_parser,
        default_value_t = false,
        help = "Do all work in memory, faster but uses more memory"
    )]
    pub in_memory: bool,

    #[clap(
        long,
        value_parser,
        default_value_t = 1,
        help = "Number of threads to use, a value of 0 uses all available threads"
    )]
    pub threads: usize,

    #[clap(
        long,
        value_parser,
        help = "Save sampled quality scores, read lengths, and insert sizes to files"
    )]
    pub save_intermediates: Option<String>,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    Generate(GenerateCommand),
    Simulate(SimulateCommand),
}

#[derive(Debug, Parser)]
#[clap(version, about, long_about = None)]
pub struct CliArgs {
    #[clap(subcommand)]
    pub command: Command,
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
