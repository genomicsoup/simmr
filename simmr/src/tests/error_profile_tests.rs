/**
 * file: error_profile_tests.rs
 * desc: Test various error profiles.
 */
use crate::error_profiles::{ErrorProfile, MinimalShortErrorProfile, PerfectShortErrorProfile};

#[test]
fn test_zero_profile() {
    let zp = PerfectShortErrorProfile {
        read_length: 150,
        insert_size: 150,
    };

    assert!(zp.get_read_length(None) == 150);
    assert!(zp.get_insert_size(None) == 150);

    let quality = zp.simulate_phred_scores(150, None);

    assert!(quality.len() == 150);
    assert!(quality.iter().all(|q| *q == 60));
}
