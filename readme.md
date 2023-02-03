
# simmr

`simmr`, short for **sim**ulating **m**etagenomic **r**eads, can be used to simulate short and 
long reads from a set of genome assemblies.
Written in Rust, `simmr` is designed to be fairly fast and efficient, and comes with a 
simple set of options for read simulation.
Although it's primarily intended to generate samples of metagenomic reads (and some of the empirical 
distributions it uses reflects this), it should work with assemblies from any organism.

## Usage

Simulate a set of 5000 perfect paired-end (PE) reads from a single genome

```bash
$ simmr --genome my-genome.fa --num-reads 5000 --output my-reads.fq
```

or from multiple genomes

```bash
$ simmr --genome genome-0.fa --genome genome-1.fa --num-reads 5000 --output my-reads.fq
```

`simmr` provides two main options for altering how reads are simulated--error and abundance profiles.
Respectively, these options configure how errors are introduced into simulated reads 
and how to handle the abundances of specific genomes in read mixtures.

### Error profiles

Two different error profiles are included with `simmr`.

The first is the `perfect` error profile. 
For short reads, this option will simulate perfect reads, i.e., reads no errors
and perfect Phred quality scores of 60. 
These reads are essentially just k-length slices of the original genome.
Long reads are simulated the same way, however individual read lengths are simulated
using a gamma distribution.

Second is the `minimal` error profile.
This profile simulates Phred scores using a normal distribution with a mean quality
score of 30.
It also assumes and implements a uniform substitution rate, i.e., the probability
of a nucleotide (e.g., A) turning into any other (e.g., C, T, or G) is 33.33%.
Substitutions only occur if the simulated quality score is lower than one randomly
generated during the substitution phase of simulation.
For long reads, read lengths are again simulated using a gamma distribution.

Use the `--error-profile` option to select a profile (the default is `perfect-short`):

```
$ simmr --error-profile perfect-short --genome my-genome.fa --output my-reads.fq
```

#### Custom error profiles

Custom error profiles, which model quality score distributions and sequencing errors, can be 
generated from real sequencing datasets.
Two custom error profiles are included with `simmr`.
The first can be used to simulate metagenomic long reads from an ONT MiniION.
The second can be used to simulate metagenomic short reads from an Illumina NoveSeq.
See the [`simmrd`](simmrd/) docs for more info on generating your own.

### Abundance profiles

There are two different profiles for controlling the abundance of genomes in a
mixture of reads with more than a single genome.
The `uniform` abundance profile uniformly distributes read abundances across all 
genomes.

```
$ simmr \
    --abundance-profile uniform \
    --num-reads 200 \
    --genome genome-0.fa \
    --genome genome-1.fa \
    --genome genome-2.fa \
    --genome genome-3.fa \
    --output reads.fq
```

e.g., simulating 200 reads from 4 genomes using the uniform profile would produce
a read mixture where the abundance of each genome was roughly 25% and therefore 50 
reads would be simulated per genome.

The `exact` abundance profile is generates `--num-reads` for each genome instead of 
distributing that number evenly like the `uniform` profile does.
In the example above, replacing `--abundance-profile uniform` with 
`--abundance-profile exact` produces a read mixture with 200 reads per genome.

### Genome file

Finally, instead of providing a list of genomes at the command line, a single file containing assembly filepaths can be provided instead.

```
$ simmr --num-reads 5000 --genome-file my-genomes.txt --output reads.fq
```

This file comes in two variants: simple and complex.
The simple variant is just a text file with filepaths to each assembly, one per line.

The complex variant is formatted as a TSV file and allows for the specification of additional genome metadata.
This file can have the following fields/columns in any order.
The `uuid` and `abundance` fields are completely optional.

- `path`
    - Genome filepath
- `uuid`
    - A unique identifier for this genome. It can be used to tag reads and track
      their provenance (which reads came from which genomes). 
- `abundance`
    - Simulate reads at this abundance percentage. If the total abundance of all 
      reads does not add up to 100.0%, the values will be normalized.

## Development

TODO

## Help options

```bash
$ simmr --help
simmr 0.1.0
Simulate metagenomic short and long reads

USAGE:
    simmr [OPTIONS] --output <OUTPUT> <--genome <GENOME>|--genome-file <GENOME_FILE>>

OPTIONS:
        --abundance-profile <ABUNDANCE_PROFILE>
            
            Genome abundances for the set of simulated reads.
            
            <full>    generates NUM_READS for each genome
            <uniform> uniformly distributes abundances across each genome
            
             [default: uniform] [possible values: full, uniform]

        --error-profile <ERROR_PROFILE>
            
            Error profile to use for read simulation
            
            <perfect-short>  simulate perfect short reads
            <perfect-long>   simulate perfect long reads
            <minimal-short>  simulate short reads, with quality scores selected from a
                             normal distribution and a uniform rate of substitution
            <minimal-long>   simulate short reads, with quality scores selected from a
                             normal distribution and a uniform rate of substitution
             [default: perfect-short] [possible values: minimal-short, minimal-long, perfect-short,
            perfect-long, ont]

        --genome <GENOME>
            Filepath to a genome to use for simulations

        --genome-file <GENOME_FILE>
            
            File containing input genome filepaths and metadata, one filepath per line
            
            If a TSV file is provided, additional metadata columns can be used for
            read simulation and abundance calculation. <fields> are all optional.
            
                filepath  <abundance>  <uuid>

    -h, --help
            Print help information

        --insert-size <INSERT_SIZE>
            Insert size for PE reads (nt) [default: 150]

        --mean-phred-score <MEAN_PHRED_SCORE>
            Average Phred quality score to use during read quality simulation [default: 30]

        --num-reads <NUM_READS>
            Number of reads to simulate [default: 1000]

        --output <OUTPUT>
            FASTQ output containing simulated reads

        --read-header-format <READ_HEADER_FORMAT>
            
            Header format for simulated reads. This string will be attached to all read headers.
            Supports
            some string interpolation via usage of the following fields:
            
                {:genome_id:}  replaced by the genome UUID the simulated read is derived from
                {:read_id:}    replaced by the unique read ID
                {:pair:}       if simulating PE reads, replaced by the pair number 1 or 2
            
             [default: @{:read_id:}|{:genome_id:}/{:pair:}]

        --read-length <READ_LENGTH>
            Individual read length (nt) [default: 150]

    -V, --version
            Print version information

        --with-ani <WITH_ANI>
            Generate reads with an average identity of N (compared to their reference)
```

## TODOs

- [x] Account for genome size when determining abundances
- [x] Additional abundance options
- [] Finish custom distribution implementation
- [] Distribution testing and plots
- [] Add `--ani` option to generate sequences with an ANI close to the reference
- [] Parallelize simulation and distribution functions
- [] Provide an option for generating separate FASTQs for each read pair (non-interleaved outputs)
- [] Provide an option for writing metadata output to a different filepath
- [] IUPAC support
- [] Additional tests
- [] Simulate inner distances and insert sizes
- [] Respect the `--mean-phred-score` parameter
- [] Handle additional CIGAR string characters.
    - Currently handles D, H, I, M, and S
- [] Flesh out readme
- [] Streaming-ish operations (as opposed to holding everything in memory then writing to disk)
- [] Tweakable long read lengths