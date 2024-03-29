/**
 * file: mod.rs
 * desc: Profile module which implements simulated read error profiles.
 */
mod base;
mod custom;
mod exact;
mod uniform;

pub use self::base::AbundanceProfile;
pub use self::custom::CustomAbundanceProfile;
pub use self::exact::ExactAbundanceProfile;
pub use self::uniform::UniformAbundanceProfile;

#[cfg(test)]
#[path = "../tests/abundance_profile_tests.rs"]
mod abundance_profile_tests;
