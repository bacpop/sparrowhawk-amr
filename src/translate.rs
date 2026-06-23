pub const DEFAULT_BACTERIAL_TRANSLATION_TABLE: u8 = 11;

pub fn translate_cds(seq: &[u8], translation_table: u8) -> Vec<u8> {
    let mut protein = Vec::with_capacity(seq.len() / 3);
    for (idx, codon) in seq.chunks_exact(3).enumerate() {
        let Some(mut aa) = translate_codon(codon, translation_table) else {
            protein.push(b'X');
            continue;
        };
        if idx == 0 && is_start_codon(codon, translation_table) {
            aa = b'M';
        }
        if aa == b'*' {
            break;
        }
        protein.push(aa);
    }
    protein
}

fn translate_codon(codon: &[u8], translation_table: u8) -> Option<u8> {
    let codon = normalize_codon(codon)?;
    let aa = match codon {
        [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => b'F',
        [b'T', b'T', b'A']
        | [b'T', b'T', b'G']
        | [b'C', b'T', b'T']
        | [b'C', b'T', b'C']
        | [b'C', b'T', b'A']
        | [b'C', b'T', b'G'] => b'L',
        [b'T', b'C', b'T']
        | [b'T', b'C', b'C']
        | [b'T', b'C', b'A']
        | [b'T', b'C', b'G']
        | [b'A', b'G', b'T']
        | [b'A', b'G', b'C'] => b'S',
        [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => b'Y',
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] => b'*',
        [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => b'C',
        [b'T', b'G', b'A'] if translation_table == 4 => b'W',
        [b'T', b'G', b'A'] => b'*',
        [b'T', b'G', b'G'] => b'W',
        [b'C', b'C', b'T'] | [b'C', b'C', b'C'] | [b'C', b'C', b'A'] | [b'C', b'C', b'G'] => b'P',
        [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => b'H',
        [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => b'Q',
        [b'C', b'G', b'T']
        | [b'C', b'G', b'C']
        | [b'C', b'G', b'A']
        | [b'C', b'G', b'G']
        | [b'A', b'G', b'A']
        | [b'A', b'G', b'G'] => b'R',
        [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => b'I',
        [b'A', b'T', b'G'] => b'M',
        [b'A', b'C', b'T'] | [b'A', b'C', b'C'] | [b'A', b'C', b'A'] | [b'A', b'C', b'G'] => b'T',
        [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => b'N',
        [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => b'K',
        [b'G', b'T', b'T'] | [b'G', b'T', b'C'] | [b'G', b'T', b'A'] | [b'G', b'T', b'G'] => b'V',
        [b'G', b'C', b'T'] | [b'G', b'C', b'C'] | [b'G', b'C', b'A'] | [b'G', b'C', b'G'] => b'A',
        [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => b'D',
        [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => b'E',
        [b'G', b'G', b'T'] | [b'G', b'G', b'C'] | [b'G', b'G', b'A'] | [b'G', b'G', b'G'] => b'G',
        _ => return None,
    };
    Some(aa)
}

fn is_start_codon(codon: &[u8], translation_table: u8) -> bool {
    let Some(codon) = normalize_codon(codon) else {
        return false;
    };
    match translation_table {
        11 => matches!(
            codon,
            [b'A', b'T', b'G'] | [b'G', b'T', b'G'] | [b'T', b'T', b'G']
        ),
        _ => matches!(codon, [b'A', b'T', b'G']),
    }
}

fn normalize_codon(codon: &[u8]) -> Option<[u8; 3]> {
    if codon.len() != 3 {
        return None;
    }
    let mut out = [b'N'; 3];
    for (idx, &base) in codon.iter().enumerate() {
        out[idx] = match base.to_ascii_uppercase() {
            b'A' => b'A',
            b'C' => b'C',
            b'G' => b'G',
            b'T' | b'U' => b'T',
            _ => return None,
        };
    }
    Some(out)
}


// ======================================================= TEST
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_bacterial_cds_and_omits_terminal_stop() {
        assert_eq!(
            translate_cds(b"ATGGCTTGGTAA", DEFAULT_BACTERIAL_TRANSLATION_TABLE),
            b"MAW"
        );
    }

    #[test]
    fn translates_common_bacterial_start_codons_as_methionine() {
        assert_eq!(
            translate_cds(b"GTGGCTTAA", DEFAULT_BACTERIAL_TRANSLATION_TABLE),
            b"MA"
        );
    }

    #[test]
    fn ambiguous_codons_become_x() {
        assert_eq!(
            translate_cds(b"ATGNNNTAA", DEFAULT_BACTERIAL_TRANSLATION_TABLE),
            b"MX"
        );
    }

    #[test]
    fn table_four_tga_is_tryptophan() {
        assert_eq!(translate_cds(b"ATGTGATAA", 4), b"MW");
    }
}
