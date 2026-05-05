#[inline(always)]
pub fn encode_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

#[inline(always)]
fn decode_base(code: u8) -> u8 {
    match code & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

pub fn reverse_complement_code(mut code: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for _ in 0..k {
        rc = (rc << 2) | ((code & 3) ^ 3);
        code >>= 2;
    }
    rc
}

#[inline(always)]
pub fn canonical_code(code: u64, k: usize) -> u64 {
    code.min(reverse_complement_code(code, k))
}

pub fn canonical_window(seq: &[u8]) -> Option<u64> {
    if seq.is_empty() || seq.len() > 31 {
        return None;
    }
    let mut code = 0u64;
    for &base in seq {
        code = (code << 2) | u64::from(encode_base(base)?);
    }
    Some(canonical_code(code, seq.len()))
}

pub fn decode_kmer(mut code: u64, k: usize) -> Vec<u8> {
    let mut out = vec![b'A'; k];
    for idx in (0..k).rev() {
        out[idx] = decode_base((code & 3) as u8);
        code >>= 2;
    }
    out
}

pub struct DnaKmerIter<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
    valid: usize,
    fwd: u64,
    rev: u64,
    mask: u64,
}

impl<'a> DnaKmerIter<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Option<Self> {
        if k == 0 || k > 31 || seq.len() < k {
            return None;
        }
        Some(Self {
            seq,
            k,
            pos: 0,
            valid: 0,
            fwd: 0,
            rev: 0,
            mask: (1u64 << (2 * k)) - 1,
        })
    }
}

impl Iterator for DnaKmerIter<'_> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let seq_pos = self.pos;
            self.pos += 1;
            let Some(base) = encode_base(self.seq[seq_pos]).map(u64::from) else {
                self.valid = 0;
                self.fwd = 0;
                self.rev = 0;
                continue;
            };
            self.fwd = ((self.fwd << 2) | base) & self.mask;
            self.rev = (self.rev >> 2) | ((base ^ 3) << (2 * (self.k - 1)));
            self.valid += 1;
            if self.valid >= self.k {
                return Some((seq_pos + 1 - self.k, self.fwd.min(self.rev)));
            }
        }
        None
    }
}

pub fn split_window(seq: &[u8]) -> Option<u64> {
    let k = seq.len();
    if k < 3 || k > 31 || k % 2 == 0 {
        return None;
    }
    let fwd = split_oriented(seq)?;
    let mut rc = Vec::with_capacity(k);
    for &base in seq.iter().rev() {
        rc.push(decode_base(encode_base(base)? ^ 3));
    }
    Some(fwd.min(split_oriented(&rc)?))
}

fn split_oriented(seq: &[u8]) -> Option<u64> {
    let half = (seq.len() - 1) / 2;
    let mut code = 0u64;
    for &base in &seq[..half] {
        code = (code << 2) | u64::from(encode_base(base)?);
    }
    for &base in &seq[half + 1..] {
        code = (code << 2) | u64::from(encode_base(base)?);
    }
    Some(code)
}

pub struct SplitKmerIter<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
}

impl<'a> SplitKmerIter<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Option<Self> {
        if k < 3 || k > 31 || k % 2 == 0 || seq.len() < k {
            return None;
        }
        Some(Self { seq, k, pos: 0 })
    }
}

impl Iterator for SplitKmerIter<'_> {
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos + self.k <= self.seq.len() {
            let pos = self.pos;
            self.pos += 1;
            if let Some(code) = split_window(&self.seq[pos..pos + self.k]) {
                return Some((pos, code));
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_matches_reverse_complement() {
        let a = canonical_window(b"ACGTT").unwrap();
        let b = canonical_window(b"AACGT").unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn skips_ambiguous_bases() {
        let kmers: Vec<_> = DnaKmerIter::new(b"ACGNACG", 3).unwrap().collect();
        assert_eq!(kmers.len(), 2);
        assert_eq!(kmers[0].0, 0);
        assert_eq!(kmers[1].0, 4);
    }

    #[test]
    fn split_kmers_skip_middle_base() {
        let a = split_window(b"AACCC").unwrap();
        let b = split_window(b"AAGCC").unwrap();
        assert_eq!(a, b);
    }
}
