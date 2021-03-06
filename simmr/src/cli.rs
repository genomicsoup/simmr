/**
 * file: cli.rs
 * desc: CLI parsing.
 */
use clap::{ArgEnum, ArgGroup, Parser};

use crate::abundance_profiles;
use crate::error_profiles;

/**
 * HELP DESCRIPTIONS
 */

static ABUNDANCE_HELP: &'static str = "
Genome abundances for the set of simulated reads.

<full>    generates NUM_READS for each genome
<uniform> uniformly distributes abundances across each genome

";

static ERROR_HELP: &'static str = "
Error profile to use for read simulation

<perfect-short>  simulate perfect short reads
<perfect-long>   simulate perfect long reads
<minimal-short>  simulate short reads, with quality scores selected from a
                 normal distribution and a uniform rate of substitution
<minimal-long>   simulate short reads, with quality scores selected from a
                 normal distribution and a uniform rate of substitution
";

static GENOME_FILE_HELP: &'static str = "
File containing input genome filepaths and metadata, one filepath per line

If a TSV file is provided, additional metadata columns can be used for
read simulation and abundance calculation. <fields> are all optional.

    filepath  <abundance>  <uuid>
";

static READ_FORMAT_HELP: &'static str = "
Header format for simulated reads. This string will be attached to all read headers. Supports 
some string interpolation via usage of the following fields:

    {:genome_id:}  replaced by the genome UUID the simulated read is derived from
    {:read_id:}    replaced by the unique read ID
    {:pair:}       if simulating PE reads, replaced by the pair number 1 or 2

";

/**
 * STRUCTS
 */

/**
 * Describes the type of error distribution to use when simulating reads.
*/
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
pub enum ErrorProfile {
    // The minimal profile only introduces substitution/point mutation errors using a uniform dist.
    MinimalShort,
    MinimalLong,
    PerfectShort,
    PerfectLong,
    ONT,
}

/**
 * Describes the type of abundance distribution to use.
*/
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
pub enum AbundanceProfile {
    // Similar to a uniform distribution but for each genome, we generate all reads, provided
    // by the --num-reads option, instead of dividing by the number of genomes
    Full,
    /// Generates a uniform abundance of all input genomes
    Uniform,
}

#[derive(Debug, Parser)]
#[clap(version, about, long_about = None)]
#[clap(group(
    ArgGroup::new("genomes")
    .required(true)
    .args(&["genome", "genome-file"])
))]
pub struct CliArgs {
    // Genome(s)
    #[clap(
        long,
        value_parser,
        help = "Filepath to a genome to use for simulations"
    )]
    pub genome: Vec<String>,

    // File of genome paths, one per line
    #[clap(
        long,
        value_parser,
        help = GENOME_FILE_HELP
    )]
    pub genome_file: Option<String>,

    // FASTQ output
    #[clap(long, value_parser, help = "FASTQ output containing simulated reads")]
    pub output: String,

    // Number of simulated reads to generate
    #[clap(
        long,
        default_value_t = 1000,
        value_parser,
        help = "Number of reads to simulate"
    )]
    pub num_reads: usize,

    // PE read length, can/will be overridden by machine-specific
    #[clap(
        long,
        default_value_t = 150,
        value_parser,
        help = "Individual read length (nt)"
    )]
    pub read_length: u16,

    // Read insert size
    #[clap(
        long,
        default_value_t = 150,
        value_parser,
        help = "Insert size for PE reads (nt) [not implemented]"
    )]
    pub insert_size: u16,

    //
    #[clap(
        long,
        default_value_t = 30,
        value_parser,
        help = "Average Phred quality score to use during read quality simulation [not implemented]"
    )]
    pub mean_phred_score: u8,

    // ErrorProfile to use for read simulation
    #[clap(
        long,
        arg_enum,
        default_value_t = ErrorProfile::PerfectShort,
        value_parser,
        help = ERROR_HELP
    )]
    pub error_profile: ErrorProfile,

    // Abundance profile to use for read simulation
    #[clap(
        long,
        arg_enum,
        default_value_t = AbundanceProfile::Uniform,
        value_parser,
        help = ABUNDANCE_HELP
    )]
    pub abundance_profile: AbundanceProfile,

    // With
    #[clap(
        long,
        value_parser,
        help = "Generate reads with an average identity of N (compared to their reference) [not implemented]"
    )]
    pub with_ani: Option<u8>,

    // Attach the following format string to all read headers
    #[clap(
        long,
        default_value = "@{:read_id:}|{:genome_id:}/{:pair:}",
        value_parser,
        help = READ_FORMAT_HELP
    )]
    pub read_header_format: String,

    #[clap(long, value_parser, help = "Random seed")]
    pub seed: Option<u64>,
}

/**
 * FUNCTIONS
 */

/**
 * Fill out and return an ErrorProfile implementation using the user provided arguments.
 */
pub fn determine_error_profile(args: &CliArgs) -> Box<dyn error_profiles::ErrorProfile> {
    match args.error_profile {
        ErrorProfile::PerfectShort => Box::new(error_profiles::PerfectShortErrorProfile {
            read_length: args.read_length,
            insert_size: args.insert_size,
        }),
        ErrorProfile::MinimalShort => Box::new(error_profiles::MinimalShortErrorProfile {
            read_length: args.read_length,
            insert_size: args.insert_size,
        }),
        ErrorProfile::PerfectLong => Box::new(error_profiles::PerfectLongErrorProfile {}),
        ErrorProfile::MinimalLong => Box::new(error_profiles::MinimalLongErrorProfile {}),
        _ => todo!(),
    }
}

/**
 * Fill out and return an AbundanceProfile implementation using the user provided arguments.
 */
pub fn determine_abundance_profile(
    args: &CliArgs,
) -> Box<dyn abundance_profiles::AbundanceProfile> {
    match args.abundance_profile {
        AbundanceProfile::Full => Box::new(abundance_profiles::UniformAbundanceProfile {}),
        AbundanceProfile::Uniform => Box::new(abundance_profiles::UniformAbundanceProfile {}),
    }
}

pub fn parse_cli_args() -> CliArgs {
    CliArgs::parse()
}
