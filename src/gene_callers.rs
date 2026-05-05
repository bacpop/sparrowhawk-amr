use crate::fasta::{FastaRecord, read_fasta};
use anyhow::{Context, ensure};
use bio::bio_types::strand::Strand;
use orphos_core::config::{OrphosConfig, OutputFormat};
use orphos_core::engine::{OrphosAnalyzer, UntrainedOrphos};
use orphos_core::output::write_results;
use orphos_core::results::OrphosResults;
use orphos_core::sequence::encoded::EncodedSequence;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

const MIN_NT_CONTIG: usize = 96;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeneCaller {
    Orphos,
}

impl std::fmt::Display for GeneCaller {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Orphos => write!(f, "orphos"),
        }
    }
}

impl std::str::FromStr for GeneCaller {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.to_ascii_lowercase().as_str() {
            "orphos" => Ok(Self::Orphos),
            _ => anyhow::bail!("unknown gene caller: {value}"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeneCallerConfig {
    pub out_dir: PathBuf,
    pub metagenomic: bool,
    pub closed_ends: bool,
    pub mask_n_runs: bool,
    pub force_non_sd: bool,
    pub translation_table: Option<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneCallerOutput {
    pub caller: GeneCaller,
    pub cds_fasta: PathBuf,
    pub protein_fasta: PathBuf,
    pub gff: PathBuf,
}

pub fn run_gene_caller(
    assembly: &Path,
    sample_name: &str,
    config: &GeneCallerConfig,
) -> anyhow::Result<GeneCallerOutput> {
    run_orphos_gene_caller(assembly, sample_name, config)
}

pub fn run_orphos_gene_caller(
    assembly: &Path,
    sample_name: &str,
    config: &GeneCallerConfig,
) -> anyhow::Result<GeneCallerOutput> {
    fs::create_dir_all(&config.out_dir)
        .with_context(|| format!("create {}", config.out_dir.display()))?;
    let cds = config.out_dir.join(format!("{sample_name}.cds.fna"));
    let proteins = config.out_dir.join(format!("{sample_name}.faa"));
    let gff = config.out_dir.join(format!("{sample_name}.gff"));

    let records = read_fasta(assembly)?;
    let results = call_orphos(&records, config)?;
    write_gff(&gff, &results)?;
    write_cds_fasta(&cds, &records, &results)?;
    fs::write(&proteins, b"")
        .with_context(|| format!("write protein placeholder {}", proteins.display()))?;

    Ok(GeneCallerOutput {
        caller: GeneCaller::Orphos,
        cds_fasta: cds,
        protein_fasta: proteins,
        gff,
    })
}

fn call_orphos(
    records: &[FastaRecord],
    config: &GeneCallerConfig,
) -> anyhow::Result<Vec<OrphosResults>> {
    let orphos_config = OrphosConfig {
        metagenomic: config.metagenomic,
        closed_ends: config.closed_ends,
        mask_n_runs: config.mask_n_runs,
        force_non_sd: config.force_non_sd,
        quiet: true,
        output_format: OutputFormat::Gff,
        translation_table: config.translation_table,
        num_threads: None,
    };
    let mut analyzer = OrphosAnalyzer::new(orphos_config.clone());

    if config.metagenomic {
        let mut results = Vec::new();
        for record in records
            .iter()
            .filter(|record| record.seq.len() >= MIN_NT_CONTIG)
        {
            results.push(
                analyzer
                    .analyze_sequence_bytes(&record.seq, record.id.clone(), description(record))
                    .map_err(|err| anyhow::anyhow!("Orphos analysis failed: {err}"))?,
            );
        }
        return Ok(results);
    }

    let mut training_seq = Vec::new();
    for record in records
        .iter()
        .filter(|record| record.seq.len() >= MIN_NT_CONTIG)
    {
        if !training_seq.is_empty() {
            training_seq.extend_from_slice(b"TTAATTAATTAA");
        }
        training_seq.extend_from_slice(&record.seq);
    }
    ensure!(
        !training_seq.is_empty(),
        "no contigs at least {MIN_NT_CONTIG} bp for Orphos"
    );

    let encoded_training = if config.mask_n_runs {
        EncodedSequence::with_masking(&training_seq)
    } else {
        EncodedSequence::without_masking(&training_seq)
    };
    let mut untrained = UntrainedOrphos::with_config(orphos_config)
        .map_err(|err| anyhow::anyhow!("Orphos configuration failed: {err}"))?;
    let training = untrained
        .train_single_genome(&encoded_training)
        .map_err(|err| anyhow::anyhow!("Orphos training failed: {err}"))?
        .into_training();

    let mut results = Vec::new();
    for record in records
        .iter()
        .filter(|record| record.seq.len() >= MIN_NT_CONTIG)
    {
        results.push(
            analyzer
                .analyze_sequence_bytes_with_training(
                    &record.seq,
                    record.id.clone(),
                    description(record),
                    training.clone(),
                )
                .map_err(|err| anyhow::anyhow!("Orphos analysis failed: {err}"))?,
        );
    }
    Ok(results)
}

fn description(record: &FastaRecord) -> Option<String> {
    record
        .description
        .strip_prefix(&record.id)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(ToString::to_string)
}

fn write_gff(path: &Path, results: &[OrphosResults]) -> anyhow::Result<()> {
    let mut out = Vec::new();
    for result in results {
        write_results(&mut out, result, OutputFormat::Gff)
            .map_err(|err| anyhow::anyhow!("write Orphos GFF: {err}"))?;
    }
    fs::write(path, out).with_context(|| format!("write {}", path.display()))
}

fn write_cds_fasta(
    path: &Path,
    records: &[FastaRecord],
    results: &[OrphosResults],
) -> anyhow::Result<()> {
    let sequence_by_id: HashMap<&str, &[u8]> = records
        .iter()
        .map(|record| (record.id.as_str(), record.seq.as_slice()))
        .collect();
    let mut out = Vec::<u8>::new();
    for result in results {
        let Some(seq) = sequence_by_id
            .get(result.sequence_info.header.as_str())
            .copied()
        else {
            continue;
        };
        for (idx, gene) in result.genes.iter().enumerate() {
            if let Some(cds) = extract_cds(
                seq,
                gene.coordinates.begin,
                gene.coordinates.end,
                gene.coordinates.strand,
            ) {
                let strand = strand_symbol(gene.coordinates.strand);
                writeln!(
                    out,
                    ">{}_{}_{}..{}_{}",
                    result.sequence_info.header,
                    idx + 1,
                    gene.coordinates.begin,
                    gene.coordinates.end,
                    strand
                )?;
                write_wrapped_fasta_seq(&mut out, &cds)?;
            }
        }
    }
    fs::write(path, out).with_context(|| format!("write {}", path.display()))
}

fn extract_cds(seq: &[u8], begin: usize, end: usize, strand: Strand) -> Option<Vec<u8>> {
    if begin == 0 || end < begin || end > seq.len() {
        return None;
    }
    let mut cds = seq[begin - 1..end].to_vec();
    if matches!(strand, Strand::Reverse) {
        reverse_complement_in_place(&mut cds);
    }
    Some(cds)
}

fn reverse_complement_in_place(seq: &mut [u8]) {
    seq.reverse();
    for base in seq {
        *base = match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            other => other,
        };
    }
}

fn strand_symbol(strand: Strand) -> char {
    match strand {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        Strand::Unknown => '.',
    }
}

fn write_wrapped_fasta_seq(out: &mut Vec<u8>, seq: &[u8]) -> anyhow::Result<()> {
    for chunk in seq.chunks(80) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_forward_cds_with_one_based_inclusive_coordinates() {
        let cds = extract_cds(b"AACCGGTT", 2, 5, Strand::Forward).unwrap();
        assert_eq!(cds, b"ACCG");
    }

    #[test]
    fn extracts_reverse_cds_as_reverse_complement() {
        let cds = extract_cds(b"AACCGGTT", 2, 5, Strand::Reverse).unwrap();
        assert_eq!(cds, b"CGGT");
    }

    #[test]
    fn rejects_out_of_bounds_coordinates() {
        assert!(extract_cds(b"AACCGGTT", 0, 5, Strand::Forward).is_none());
        assert!(extract_cds(b"AACCGGTT", 5, 2, Strand::Forward).is_none());
        assert!(extract_cds(b"AACCGGTT", 1, 99, Strand::Forward).is_none());
    }
}
