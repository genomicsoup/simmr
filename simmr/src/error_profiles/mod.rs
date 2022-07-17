/**
 * file: mod.rs
 * desc: Profile module which implements simulated read error profiles.
 */
mod base;
mod custom_long;
mod minimal_long;
mod minimal_short;
mod perfect_long;
mod perfect_short;

pub use self::base::ErrorProfile;
pub use self::minimal_long::MinimalLongErrorProfile;
pub use self::minimal_short::MinimalShortErrorProfile;
pub use self::perfect_long::PerfectLongErrorProfile;
pub use self::perfect_short::PerfectShortErrorProfile;

#[cfg(test)]
#[path = "../tests/error_profile_tests.rs"]
mod error_profile_tests;
