use crate::amrfinder_db::load_amrfinder_references;
use crate::detect::DetectionResult;
use crate::fasta::read_fasta;
use crate::index::{AmrIndex, KmerAssignment};
use crate::kmer::{DnaKmerIter, SplitKmerIter, decode_kmer};
use anyhow::{Context, ensure};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct DebugMissesConfig<'a> {
    pub index: &'a AmrIndex,
    pub assembly_path: &'a Path,
    pub amrfinder_tsv: &'a Path,
    pub detector_json: &'a Path,
    pub db_dir: Option<&'a Path>,
    pub refinement_k: usize,
    pub missing_kmer_limit: usize,
}

#[derive(Debug, Serialize)]
pub struct DebugMissesReport {
    pub assembly: String,
    pub detector_json: String,
    pub amrfinder_tsv: String,
    pub index_k: usize,
    pub refinement_k: usize,
    pub amrfinder_amr_gene_count: usize,
    pub detector_gene_exact_count: usize,
    pub missed_gene_count: usize,
    pub missed: Vec<DebugMiss>,
}

#[derive(Debug, Serialize)]
pub struct DebugMiss {
    pub element_symbol: String,
    pub hierarchy_node: String,
    pub method: String,
    pub contig_id: String,
    pub start: usize,
    pub stop: usize,
    pub strand: String,
    pub target_length: Option<usize>,
    pub reference_sequence_length: Option<usize>,
    pub coverage_of_reference: Option<f64>,
    pub identity_to_reference: Option<f64>,
    pub closest_reference_accession: String,
    pub index_gene_found_by_accession: bool,
    pub index_gene_id: Option<String>,
    pub index_element_symbol: Option<String>,
    pub index_gene_symbol: Option<String>,
    pub index_allele_symbol: Option<String>,
    pub index_family: Option<String>,
    pub index_hierarchy_node: Option<String>,
    pub index_gene_length: Option<usize>,
    pub interval_length: usize,
    pub interval_distinct_k31: usize,
    pub gene_diagnostic_k31_total: usize,
    pub gene_diagnostic_k31_matched: usize,
    pub gene_diagnostic_k31_missing: usize,
    pub gene_diagnostic_k31_fraction: f64,
    pub detector_gene_threshold_met: bool,
    pub detector_seed_threshold_met: bool,
    pub family_diagnostic_k31_total: usize,
    pub family_diagnostic_k31_matched: usize,
    pub family_diagnostic_k31_fraction: f64,
    pub detector_family_threshold_met: bool,
    pub lowk_refinement_target_total: usize,
    pub lowk_refinement_matched: usize,
    pub lowk_refinement_fraction: f64,
    pub split_refinement_target_total: usize,
    pub split_refinement_matched: usize,
    pub split_refinement_fraction: f64,
    pub reference_unique_k31_total: Option<usize>,
    pub reference_unique_k31_matched: Option<usize>,
    pub reference_unique_k31_fraction: Option<f64>,
    pub missing_diagnostic_k31_examples: Vec<String>,
}

#[derive(Debug, Clone)]
struct AmrfinderRow {
    element_symbol: String,
    hierarchy_node: String,
    method: String,
    contig_id: String,
    start: usize,
    stop: usize,
    strand: String,
    target_length: Option<usize>,
    reference_sequence_length: Option<usize>,
    coverage_of_reference: Option<f64>,
    identity_to_reference: Option<f64>,
    closest_reference_accession: String,
}

pub fn debug_amrfinder_misses(config: DebugMissesConfig<'_>) -> anyhow::Result<DebugMissesReport> {
    ensure!(
        config.refinement_k > 0 && config.refinement_k <= config.index.k,
        "refinement k must be between 1 and index k"
    );

    let assembly = read_fasta(config.assembly_path)?;
    let contigs: HashMap<String, Vec<u8>> = assembly
        .into_iter()
        .map(|record| (record.id, record.seq))
        .collect();
    let amrfinder_rows = parse_amrfinder_rows(config.amrfinder_tsv)?;
    let detector = parse_detector_json(config.detector_json)?;
    let detector_gene_symbols = detector_gene_symbols(&detector);
    let reference_kmers = if let Some(db_dir) = config.db_dir {
        Some(reference_kmers_by_accession(db_dir, config.index.k)?)
    } else {
        None
    };
    let gene_by_accession: HashMap<&str, usize> = config
        .index
        .genes
        .iter()
        .enumerate()
        .map(|(idx, gene)| (gene.protein_accession.as_str(), idx))
        .collect();

    let mut missed = Vec::new();
    for row in &amrfinder_rows {
        if detector_gene_symbols.contains(&row.element_symbol) {
            continue;
        }
        missed.push(debug_row(
            row,
            config.index,
            &gene_by_accession,
            &contigs,
            reference_kmers.as_ref(),
            config.refinement_k,
            config.missing_kmer_limit,
        )?);
    }

    Ok(DebugMissesReport {
        assembly: config.assembly_path.display().to_string(),
        detector_json: config.detector_json.display().to_string(),
        amrfinder_tsv: config.amrfinder_tsv.display().to_string(),
        index_k: config.index.k,
        refinement_k: config.refinement_k,
        amrfinder_amr_gene_count: amrfinder_rows.len(),
        detector_gene_exact_count: detector_gene_symbols.len(),
        missed_gene_count: missed.len(),
        missed,
    })
}

fn debug_row(
    row: &AmrfinderRow,
    index: &AmrIndex,
    gene_by_accession: &HashMap<&str, usize>,
    contigs: &HashMap<String, Vec<u8>>,
    reference_kmers: Option<&HashMap<String, HashSet<u64>>>,
    refinement_k: usize,
    missing_kmer_limit: usize,
) -> anyhow::Result<DebugMiss> {
    let interval = extract_interval(contigs, row)?;
    let interval_kmers = distinct_kmers(&interval, index.k);
    let gene_id = gene_by_accession
        .get(row.closest_reference_accession.as_str())
        .copied();

    let (
        index_gene_id,
        index_element_symbol,
        index_gene_symbol,
        index_allele_symbol,
        index_family,
        index_hierarchy_node,
        index_gene_length,
        diagnostic_total,
        diagnostic_matched,
        family_diagnostic_total,
        family_diagnostic_matched,
        missing_diagnostic,
        lowk_total,
        lowk_matched,
        split_total,
        split_matched,
    ) = if let Some(gene_id) = gene_id {
        let gene = &index.genes[gene_id];
        let diagnostic: HashSet<u64> = index.gene_specific_kmers(gene_id).iter().copied().collect();
        let matched: HashSet<u64> = diagnostic.intersection(&interval_kmers).copied().collect();
        let missing: Vec<u64> = diagnostic.difference(&interval_kmers).copied().collect();
        let family_diagnostic = family_specific_kmers(index, &gene.family);
        let family_matched = family_diagnostic.intersection(&interval_kmers).count();
        let lowk_target = refinement_target(&missing, index.k, refinement_k, false);
        let split_target = refinement_target(&missing, index.k, refinement_k, true);
        let interval_lowk = distinct_lowk(&interval, refinement_k);
        let interval_split = distinct_split(&interval, refinement_k);
        (
            Some(gene.id.clone()),
            Some(gene.element_symbol.clone()),
            Some(gene.gene_symbol.clone()),
            Some(gene.allele_symbol.clone()),
            Some(gene.family.clone()),
            Some(gene.hierarchy_node.clone()),
            Some(gene.length),
            diagnostic.len(),
            matched.len(),
            family_diagnostic.len(),
            family_matched,
            missing,
            lowk_target.len(),
            lowk_target.intersection(&interval_lowk).count(),
            split_target.len(),
            split_target.intersection(&interval_split).count(),
        )
    } else {
        (
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            0,
            0,
            0,
            0,
            Vec::new(),
            0,
            0,
            0,
            0,
        )
    };

    let diagnostic_fraction = fraction(diagnostic_matched, diagnostic_total);
    let family_fraction = fraction(family_diagnostic_matched, family_diagnostic_total);
    let reference_counts = reference_kmers
        .and_then(|by_accession| by_accession.get(&row.closest_reference_accession))
        .map(|kmers| {
            let matched = kmers.intersection(&interval_kmers).count();
            (kmers.len(), matched, fraction(matched, kmers.len()))
        });

    Ok(DebugMiss {
        element_symbol: row.element_symbol.clone(),
        hierarchy_node: row.hierarchy_node.clone(),
        method: row.method.clone(),
        contig_id: row.contig_id.clone(),
        start: row.start,
        stop: row.stop,
        strand: row.strand.clone(),
        target_length: row.target_length,
        reference_sequence_length: row.reference_sequence_length,
        coverage_of_reference: row.coverage_of_reference,
        identity_to_reference: row.identity_to_reference,
        closest_reference_accession: row.closest_reference_accession.clone(),
        index_gene_found_by_accession: gene_id.is_some(),
        index_gene_id,
        index_element_symbol,
        index_gene_symbol,
        index_allele_symbol,
        index_family,
        index_hierarchy_node,
        index_gene_length,
        interval_length: interval.len(),
        interval_distinct_k31: interval_kmers.len(),
        gene_diagnostic_k31_total: diagnostic_total,
        gene_diagnostic_k31_matched: diagnostic_matched,
        gene_diagnostic_k31_missing: diagnostic_total.saturating_sub(diagnostic_matched),
        gene_diagnostic_k31_fraction: diagnostic_fraction,
        detector_gene_threshold_met: diagnostic_fraction >= 0.10,
        detector_seed_threshold_met: diagnostic_matched >= 3 || diagnostic_fraction >= 0.01,
        family_diagnostic_k31_total: family_diagnostic_total,
        family_diagnostic_k31_matched: family_diagnostic_matched,
        family_diagnostic_k31_fraction: family_fraction,
        detector_family_threshold_met: family_fraction >= 0.10,
        lowk_refinement_target_total: lowk_total,
        lowk_refinement_matched: lowk_matched,
        lowk_refinement_fraction: fraction(lowk_matched, lowk_total),
        split_refinement_target_total: split_total,
        split_refinement_matched: split_matched,
        split_refinement_fraction: fraction(split_matched, split_total),
        reference_unique_k31_total: reference_counts.map(|counts| counts.0),
        reference_unique_k31_matched: reference_counts.map(|counts| counts.1),
        reference_unique_k31_fraction: reference_counts.map(|counts| counts.2),
        missing_diagnostic_k31_examples: missing_diagnostic
            .into_iter()
            .take(missing_kmer_limit)
            .map(|kmer| String::from_utf8_lossy(&decode_kmer(kmer, index.k)).to_string())
            .collect(),
    })
}

fn parse_detector_json(path: &Path) -> anyhow::Result<DetectionResult> {
    let text = fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
    serde_json::from_str(&text).with_context(|| format!("parse {}", path.display()))
}

fn detector_gene_symbols(result: &DetectionResult) -> HashSet<String> {
    result
        .hits
        .iter()
        .filter(|hit| hit.call_type == "gene")
        .filter_map(|hit| hit.element_symbol.clone())
        .filter(|symbol| !symbol.is_empty())
        .collect()
}

fn parse_amrfinder_rows(path: &Path) -> anyhow::Result<Vec<AmrfinderRow>> {
    let text = fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let header = lines.next().context("AMRFinderPlus TSV is empty")?;
    let columns: HashMap<&str, usize> = header
        .split('\t')
        .enumerate()
        .map(|(i, c)| (c, i))
        .collect();
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let type_value = tsv_field(&fields, &columns, "Type");
        let subtype_value = tsv_field(&fields, &columns, "Subtype");
        if type_value != "AMR" || subtype_value != "AMR" {
            continue;
        }
        rows.push(AmrfinderRow {
            element_symbol: tsv_field(&fields, &columns, "Element symbol").to_string(),
            hierarchy_node: tsv_field(&fields, &columns, "Hierarchy node").to_string(),
            method: tsv_field(&fields, &columns, "Method").to_string(),
            contig_id: tsv_field(&fields, &columns, "Contig id").to_string(),
            start: parse_usize(tsv_field(&fields, &columns, "Start")),
            stop: parse_usize(tsv_field(&fields, &columns, "Stop")),
            strand: tsv_field(&fields, &columns, "Strand").to_string(),
            target_length: parse_optional_usize(tsv_field(&fields, &columns, "Target length")),
            reference_sequence_length: parse_optional_usize(tsv_field(
                &fields,
                &columns,
                "Reference sequence length",
            )),
            coverage_of_reference: parse_optional_f64(tsv_field(
                &fields,
                &columns,
                "% Coverage of reference",
            )),
            identity_to_reference: parse_optional_f64(tsv_field(
                &fields,
                &columns,
                "% Identity to reference",
            )),
            closest_reference_accession: tsv_field(
                &fields,
                &columns,
                "Closest reference accession",
            )
            .to_string(),
        });
    }
    Ok(rows)
}

fn tsv_field<'a>(fields: &'a [&str], columns: &HashMap<&str, usize>, column: &str) -> &'a str {
    columns
        .get(column)
        .and_then(|idx| fields.get(*idx))
        .copied()
        .unwrap_or("")
        .trim()
}

fn parse_usize(value: &str) -> usize {
    value.parse().unwrap_or(0)
}

fn parse_optional_usize(value: &str) -> Option<usize> {
    if value.is_empty() || value == "NA" {
        None
    } else {
        value.parse().ok()
    }
}

fn parse_optional_f64(value: &str) -> Option<f64> {
    if value.is_empty() || value == "NA" {
        None
    } else {
        value.parse().ok()
    }
}

fn extract_interval(
    contigs: &HashMap<String, Vec<u8>>,
    row: &AmrfinderRow,
) -> anyhow::Result<Vec<u8>> {
    let contig = contigs
        .get(&row.contig_id)
        .with_context(|| format!("contig {} not found in assembly FASTA", row.contig_id))?;
    let start = row.start.min(row.stop);
    let stop = row.start.max(row.stop);
    ensure!(
        start >= 1 && stop <= contig.len(),
        "AMRFinderPlus interval {}:{}-{} outside contig length {}",
        row.contig_id,
        row.start,
        row.stop,
        contig.len()
    );
    Ok(contig[start - 1..stop].to_vec())
}

fn distinct_kmers(seq: &[u8], k: usize) -> HashSet<u64> {
    DnaKmerIter::new(seq, k)
        .map(|iter| iter.map(|(_, kmer)| kmer).collect())
        .unwrap_or_default()
}

fn distinct_lowk(seq: &[u8], k: usize) -> HashSet<u64> {
    DnaKmerIter::new(seq, k)
        .map(|iter| iter.map(|(_, kmer)| kmer).collect())
        .unwrap_or_default()
}

fn distinct_split(seq: &[u8], k: usize) -> HashSet<u64> {
    SplitKmerIter::new(seq, k)
        .map(|iter| iter.map(|(_, kmer)| kmer).collect())
        .unwrap_or_default()
}

fn refinement_target(
    missing_kmers: &[u64],
    index_k: usize,
    refinement_k: usize,
    split: bool,
) -> HashSet<u64> {
    let mut target = HashSet::new();
    for &kmer in missing_kmers {
        let decoded = decode_kmer(kmer, index_k);
        if split {
            if let Some(iter) = SplitKmerIter::new(&decoded, refinement_k) {
                target.extend(iter.map(|(_, code)| code));
            }
        } else if let Some(iter) = DnaKmerIter::new(&decoded, refinement_k) {
            target.extend(iter.map(|(_, code)| code));
        }
    }
    target
}

fn family_specific_kmers(index: &AmrIndex, family: &str) -> HashSet<u64> {
    let Some(family_id) = index
        .families
        .iter()
        .position(|indexed_family| indexed_family == family)
    else {
        return HashSet::new();
    };
    index
        .kmer_codes
        .iter()
        .copied()
        .filter(|&kmer| index.lookup(kmer) == Some(KmerAssignment::Family(family_id)))
        .collect()
}

fn reference_kmers_by_accession(
    db_dir: &Path,
    k: usize,
) -> anyhow::Result<HashMap<String, HashSet<u64>>> {
    let mut by_accession = HashMap::<String, HashSet<u64>>::new();
    for reference in load_amrfinder_references(db_dir)? {
        by_accession
            .entry(reference.protein_accession)
            .or_default()
            .extend(distinct_kmers(&reference.seq, k));
    }
    Ok(by_accession)
}

fn fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}
