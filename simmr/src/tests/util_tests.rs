/**
 * file: util_tests.rs
 * desc: Utility function tests.
 */
use crate::util;

#[test]
fn test_complement_0() {
    let s1 = "aacctg".as_bytes();

    assert!(
        std::str::from_utf8(
            &s1.into_iter()
                .map(|n| util::complement(*n))
                .collect::<Vec<u8>>()
        )
        .unwrap()
            == "ttggac"
    );
}

#[test]
fn test_complement_1() {
    let s1 = "TAGCNNNN".as_bytes();

    assert!(
        std::str::from_utf8(
            &s1.into_iter()
                .map(|n| util::complement(*n))
                .collect::<Vec<u8>>()
        )
        .unwrap()
            == "ATCGNNNN"
    );
}

#[test]
fn test_complement_2() {
    let s1 = "CaTTagG".as_bytes();

    assert!(
        std::str::from_utf8(
            &s1.into_iter()
                .map(|n| util::complement(*n))
                .collect::<Vec<u8>>()
        )
        .unwrap()
            == "GtAAtcC"
    );
}

#[test]
fn test_generate_id_using_seed() {
    let bytes: [u8; 16] = [
        0xDE, 0xAD, 0xBE, 0xEF, 0xBA, 0xAD, 0xF0, 0x0D, 0xDE, 0xAD, 0xBE, 0xEF, 0xBA, 0xAD, 0xF0,
        0x0D,
    ];

    let uuid1 = util::generate_id(Some(&bytes));
    let uuid2 = util::generate_id(Some(&bytes));
    let uuid3 = util::generate_id(None);

    assert!(uuid1 == uuid2);
    assert!(uuid1 != uuid3);
    assert!(uuid2 != uuid3);
}

#[test]
fn test_encode_quality_score() {
    assert!(util::encode_quality_score(0) == b'!');
    assert!(util::encode_quality_score(1) == b'"');
    assert!(util::encode_quality_score(10) == b'+');
    assert!(util::encode_quality_score(41) == b'J');
}

#[test]
fn test_encode_quality_scores() {
    let scores = vec![0, 1, 10, 41];

    assert!(util::encode_quality_scores(&scores) == b"!\"+J");
}

#[test]
fn test_probability_to_phred_conversions() {
    assert!(util::convert_probability_to_phred(0.1) == 10);
    assert!(util::convert_probability_to_phred(0.001) == 30);
    assert!(util::convert_probability_to_phred(0.000001) == 60);
}

#[test]
fn test_phred_to_probability_conversions() {
    assert!(util::convert_phred_to_probability(10) == 0.1);
    assert!(util::convert_phred_to_probability(30) == 0.001);
    assert!(util::convert_phred_to_probability(60) == 0.000001);
}

#[test]
fn test_phred_to_accuracy_conversion() {
    assert!(util::convert_phred_to_accuracy(10) == 0.9);
    assert!(util::convert_phred_to_accuracy(30) == 0.999);
    assert!(util::convert_phred_to_accuracy(60) == 0.999999);
}

#[test]
fn test_accuracy_to_phred_conversion() {
    assert!(util::convert_accuracy_to_phred(0.9) == 10);
    assert!(util::convert_accuracy_to_phred(0.999) == 30);
    assert!(util::convert_accuracy_to_phred(0.999999) == 60);
}
