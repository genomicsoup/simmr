
# simmrd

`simmrd`, short for **sim**ulating **m**etagenomic **r**ead **d**istributions, can be used to
generate distributions of quality scores, sequencing errors (substitutions, indels), and 
read lengths from real sequencing data.

## Usage

`simmrd` requires two arguments: a SAM file containing alignments to one or more reference 
genomes (`--sam-file`) and an output file to store the serialized sequencing model (`--output`).

```bash
$ simmrd --sam-file alignments.sam --output error-model.bin
```

## How it works

### Quality scores

### Sequencing errors


## Development

TODO

## Help options

```bash
$ simmrd --help
simmrd 0.1.0
Generate distributions for simulating metagenomic short and long reads

USAGE:
    simmrd [OPTIONS] --sam-file <SAM_FILE> --output <OUTPUT>

OPTIONS:
        --bin-size <BIN_SIZE>
            Quality score bin size used for kernel density estimation [default: 5]

    -h, --help
            Print help information

        --k <K>
            Kmer length, must be between 3 and 10 [default: 7]

        --mapq-threshold <MAPQ_THRESHOLD>
            MAPQ threshold, alignments below the threshold will not be used

        --max-alignments <MAX_ALIGNMENTS>
            Use a maximum of N alignments for distribution modeling

        --max-alt-kmers <MAX_ALT_KMERS>
            Use a maximum of N alternately sequenced kmers for distribution modeling [default: 20]

        --output <OUTPUT>
            Output file

        --sam-file <SAM_FILE>
            SAM file

    -V, --version
            Print version information

```

## TODOs

- [] Handle BAM files
- [] Performance optimizations