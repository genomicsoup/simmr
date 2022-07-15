/**
 * file: util.rs
 * desc: Misc. utility functions.
 */

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
}
