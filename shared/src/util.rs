/**
 * file: util.rs
 * desc: Misc. utility functions.
 */
use num_traits::Float;
use std::iter::Sum;

/**
 * Sequence related functions
 */

/**
 * Generate the complement of the given nucleotide.
 * TODO: IUPAC support
 */
pub fn complement(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'a' => b't',
        b'T' => b'A',
        b't' => b'a',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        x => x,
    }
}

/**
 * Given a sequence of nucleotides, generate the reverse complement.
 */
pub fn reverse_complement(nucs: &[u8]) -> Vec<u8> {
    nucs.iter()
        .rev()
        .map(|n| complement(*n))
        .collect::<Vec<u8>>()
}

/**
 * Convert a byte array into a string.
 *
 * args
 *  bs: byte array
 *
 * returns
 *  a string representation of the given bytes or an empty string if the conversion failed
 */
pub fn bytes_to_string(bs: &[u8]) -> String {
    std::str::from_utf8(bs).unwrap_or("").to_string()
}

/**
 * maths
 */

pub fn mean2<T, I: Iterator<Item = T>>(iter: I) -> f64
where
    T: Into<f64> + Sum<T>,
{
    let mut length = 0;
    let sum = iter
        .map(|t| {
            length += 1;
            t
        })
        .sum::<T>();

    if length == 0 {
        return 0.0;
    }

    sum.into() / length as f64
}

pub fn mean<T>(vs: &[T]) -> T
where
    T: Float,
{
    //vs.iter().map(|t| *t).sum::<T>() / num_traits::cast(vs.len()).unwrap()
    vs.iter().fold(T::zero(), |ac: T, v| ac + *v) / num_traits::cast(vs.len()).unwrap()
}

pub fn variance<T>(vs: &[T]) -> T
where
    T: Float,
{
    let avg = mean(vs);

    vs.iter()
        .fold(T::zero(), |ac: T, v| ac + (*v - avg) * (*v - avg))
        / num_traits::cast(vs.len()).unwrap()
}

pub fn std_deviation<T>(vs: &[T]) -> T
where
    T: Float,
{
    variance(vs).sqrt()
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_complement_0() {
        let s1 = "aacctg".as_bytes();

        assert!(
            std::str::from_utf8(&s1.into_iter().map(|n| complement(*n)).collect::<Vec<u8>>())
                .unwrap()
                == "ttggac"
        );
    }

    #[test]
    fn test_complement_1() {
        let s1 = "TAGCNNNN".as_bytes();

        assert!(
            std::str::from_utf8(&s1.into_iter().map(|n| complement(*n)).collect::<Vec<u8>>())
                .unwrap()
                == "ATCGNNNN"
        );
    }

    #[test]
    fn test_complement_2() {
        let s1 = "CaTTagG".as_bytes();

        assert!(
            std::str::from_utf8(&s1.into_iter().map(|n| complement(*n)).collect::<Vec<u8>>())
                .unwrap()
                == "GtAAtcC"
        );
    }

    #[test]
    fn test_bytes_to_string() {
        assert!(bytes_to_string(&vec![b'f', b'o', b'o']) == "foo".to_string());
    }

    #[test]
    fn test_mean_1() {
        assert!(mean(&vec![1.0, 2.0, 3.0, 4.0, 5.0]) == 3.0);
    }

    #[test]
    fn test_std_dev_1() {
        assert!((std_deviation(&vec![1.0, 2.0, 3.0, 4.0, 5.0]) * 100.0).round() / 100.0 == 1.41);
    }
}
