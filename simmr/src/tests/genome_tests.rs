/**
 * file: genome_tests.rs
 * desc: Genome struct and fasta parsing tests.
 */
use crate::genome;

#[test]
fn test_genome_from_fasta() {
    let genome = genome::Genome::from_fasta("src/tests/data/sample.fna").unwrap();

    assert!(genome.filepath.to_str().unwrap() == "src/tests/data/sample.fna");
    assert!(genome.num_seqs == 2);
    assert!(genome.size == 320);

    assert!(std::str::from_utf8(&genome.sequence[0].id).unwrap() == "header1");
    assert!(genome.sequence[0].size == 160);

    assert!(std::str::from_utf8(&genome.sequence[1].id).unwrap() == "header2");
    assert!(genome.sequence[1].size == 160);
}
