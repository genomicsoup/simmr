/**
 * file: abundance_profile_tests.rs
 * desc: Test various abundance profiles.
 */
use crate::abundance_profiles::{AbundanceProfile, UniformAbundanceProfile};

#[test]
fn test_uniform_profile() {
    let up = UniformAbundanceProfile { size_aware: false };
    let abundances = up.determine_abundances(100, 5);

    assert!(abundances.len() == 5);

    assert!(abundances[0].0 == 20);
    assert!(abundances[0].1 == 20.0);

    assert!(abundances[1].0 == 20);
    assert!(abundances[1].1 == 20.0);

    assert!(abundances[2].0 == 20);
    assert!(abundances[2].1 == 20.0);

    assert!(abundances[3].0 == 20);
    assert!(abundances[3].1 == 20.0);

    assert!(abundances[4].0 == 20);
    assert!(abundances[4].1 == 20.0);
}
