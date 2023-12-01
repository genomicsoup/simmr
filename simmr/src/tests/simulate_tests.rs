/**
 * file: simulate_tests.rs
 * desc: Test functions used for simulating reads.
 */
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

use crate::error_profiles;
use crate::genome;
use crate::simulate;

#[test]
#[ignore = "updated the RNG and have to redo ground truth"]
fn test_simulate_pe_read_using_zero_error_profile() {
    // The no error profile
    let error_profile = error_profiles::PerfectShortErrorProfile {
        read_length: 20,
        insert_size: 20,
    };
    let genome = genome::Genome::from_fasta(
        "src/tests/data/GCF_000005845.2_ASM584v2_genomic.partial.fna",
        false,
    )
    .unwrap();

    // With this random seed, the generated seed should be 97132697663989775522 and the
    // read positions should be:
    // fwd start: 38, fwd end: 58, rev start: 98, rev end: 78
    let mut rng = StdRng::seed_from_u64(42);
    let seed: u64 = rng.gen();

    let pe_read =
        simulate::simulate_pe_read(&genome.sequence[0], &error_profile, Some(seed)).unwrap();

    // forward read, TGTGGATTAAAAAAAGAGTG
    assert!(std::str::from_utf8(&pe_read.forward.sequence).unwrap() == "TGTGGATTAAAAAAAGAGTG");
    // reverse read, TGGTTACCTGCCGTGAGTAA
    // reverse read complement, ACCAATGGACGGCACTCATT
    // reverse read reverse complement, TTACTCACGGCAGGTAACCA
    //assert!(std::str::from_utf8(&pe_read.reverse.sequence).unwrap() == "ACCAATGGACGGCACTCATT");
    assert!(
        std::str::from_utf8(&pe_read.reverse.as_ref().unwrap().sequence).unwrap()
            == "TTACTCACGGCAGGTAACCA"
    );

    // "perfect" quality scores
    assert!(pe_read.forward.quality.iter().all(|q| *q == 60));
    assert!(pe_read
        .reverse
        .as_ref()
        .unwrap()
        .quality
        .iter()
        .all(|q| *q == 60));
    //assert!(std::str::from_utf8(&pe_read.id)
    //    .unwrap()
    //    .contains("R0:seq_id="));
}

#[test]
#[ignore = "updated the RNG and have to redo ground truth"]
fn test_simulate_reads_using_zero_error_profile() {
    // The no error profile
    let error_profile = error_profiles::PerfectShortErrorProfile {
        read_length: 20,
        insert_size: 20,
    };
    let genome = genome::Genome::from_fasta(
        "src/tests/data/GCF_000005845.2_ASM584v2_genomic.partial.fna",
        false,
    )
    .unwrap();

    // With this random seed, the gen seed used by simulate_reads should be 6335..6202 and the
    // read positions should be:
    // fwd start: 8, fwd end: 28, rev start: 68, rev end: 48
    let reads =
        simulate::simulate_pe_reads_from_genome(1, &genome, &error_profile, Some(42)).unwrap();
    let pe_read = &reads[0];

    // forward read, ATTCTGACTGCAACGGGCAA
    assert!(std::str::from_utf8(&pe_read.forward.sequence).unwrap() == "ATTCTGACTGCAACGGGCAA");
    // reverse read, AAAAAGAGTGTCTGATAGCA
    // reverse read complement, TCTGAACTGGTTACCTGCCG
    // reverse read reverse complement, TGCTATCAGACACTCTTTTT
    //assert!(std::str::from_utf8(&pe_read.reverse.sequence).unwrap() == "TTTTTCTCACAGACTATCGT");
    assert!(
        std::str::from_utf8(&pe_read.reverse.as_ref().unwrap().sequence).unwrap()
            == "TGCTATCAGACACTCTTTTT"
    );

    // "perfect" quality scores
    assert!(pe_read.forward.quality.iter().all(|q| *q == 60));
    assert!(pe_read
        .reverse
        .as_ref()
        .unwrap()
        .quality
        .iter()
        .all(|q| *q == 60));
    //assert!(std::str::from_utf8(&pe_read.id)
    //    .unwrap()
    //    .contains("R0:seq_id="));
}
