use anyhow::{Context, ensure};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub seq: Vec<u8>,
}

pub fn read_fasta(path: &Path) -> anyhow::Result<Vec<FastaRecord>> {
    let bytes = fs::read(path).with_context(|| format!("read FASTA {}", path.display()))?;
    parse_fasta_bytes(&bytes)
}

// very simple parser, probably could do something to support compressed things like the things we've got in sketchlib or sphk-asm
pub fn parse_fasta_bytes(bytes: &[u8]) -> anyhow::Result<Vec<FastaRecord>> {
    let text = std::str::from_utf8(bytes).context("FASTA is not valid UTF-8")?;
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_description = String::new();
    let mut current_seq = Vec::new();

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push(FastaRecord {
                    id,
                    description: std::mem::take(&mut current_description),
                    seq: std::mem::take(&mut current_seq),
                });
            }
            let id = rest
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            current_id = Some(id);
            current_description = rest.to_string();
        } else {
            current_seq.extend_from_slice(trimmed.as_bytes());
        }
    }

    if let Some(id) = current_id {
        records.push(FastaRecord {
            id,
            description: current_description,
            seq: current_seq,
        });
    }

    ensure!(!records.is_empty(), "no FASTA records found");
    Ok(records)
}


// =================== TEST
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multiple_records() {
        let records = parse_fasta_bytes(b">a one\nAC\nGT\n>b\nTT\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "a");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].seq, b"TT");
    }
}
