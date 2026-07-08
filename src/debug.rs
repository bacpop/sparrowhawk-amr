use crate::amrfinder_db::{ReferenceType, load_amrfinder_references};
use crate::detect::DetectionResult;
use crate::fasta::read_fasta;
use crate::index::{AmrIndex, IndexAlphabet, ReportUnitKind, UnitId};
use crate::kmer::{DnaKmerIter, SplitKmerIter, decode_kmer};
use anyhow::{Context, ensure};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

// Helpers for debugging and analysing numbers of k-mers etc.

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

#[derive(Debug, Clone)]
pub struct TruthKmerEvidenceConfig<'a> {
    pub index: &'a AmrIndex,
    pub assembly_path: &'a Path,
    pub amrfinder_tsv: &'a Path,
    pub detector_json: &'a Path,
    pub include_types: &'a [ReferenceType],
    pub min_gene_fraction: f64,
    pub min_family_fraction: f64,
}

#[derive(Debug, Serialize)]
pub struct TruthKmerEvidenceReport {
    pub assembly: String,
    pub detector_json: String,
    pub amrfinder_tsv: String,
    pub index_k: usize,
    pub truth_count: usize,
    pub rows: Vec<TruthKmerEvidenceRow>,
}

#[derive(Debug, Serialize)]
pub struct TruthKmerEvidenceRow {
    pub element_symbol: String,
    pub hierarchy_node: String,
    pub method: String,
    pub type_name: String,
    pub subtype: String,
    pub class_name: String,
    pub subclass: String,
    pub contig_id: String,
    pub start: usize,
    pub stop: usize,
    pub strand: String,
    pub coverage_of_reference: Option<f64>,
    pub identity_to_reference: Option<f64>,
    pub closest_reference_accession: String,
    pub covered_by_detector: bool,
    pub covering_detector_units: Vec<String>,
    pub truth_supported_by_index: bool,
    pub best_index_unit: Option<String>,
    pub best_index_unit_type: Option<String>,
    pub best_index_unit_label: Option<String>,
    pub best_diagnostic_total: usize,
    pub best_diagnostic_matched: usize,
    pub best_diagnostic_missing: usize,
    pub best_diagnostic_fraction: f64,
    pub exact_diagnostic_total: usize,
    pub exact_diagnostic_matched: usize,
    pub exact_diagnostic_fraction: f64,
    pub family_diagnostic_total: usize,
    pub family_diagnostic_matched: usize,
    pub family_diagnostic_fraction: f64,
    pub interval_length: usize,
    pub interval_distinct_kmers: usize,
    pub recall_failure_category: String,
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
    type_name: String,
    subtype: String,
    class_name: String,
    subclass: String,
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

pub fn truth_kmer_evidence(
    config: TruthKmerEvidenceConfig<'_>,
) -> anyhow::Result<TruthKmerEvidenceReport> {
    ensure!(
        config.index.alphabet == IndexAlphabet::Dna,
        "AMRFinder truth k-mer evidence only supports DNA indexes"
    );

    let assembly = read_fasta(config.assembly_path)?;
    let contigs: HashMap<String, Vec<u8>> = assembly
        .into_iter()
        .map(|record| (record.id, record.seq))
        .collect();
    let include_types: HashSet<&str> = config
        .include_types
        .iter()
        .map(|value| value.as_str())
        .collect();
    let amrfinder_rows = parse_amrfinder_rows_with_types(config.amrfinder_tsv, &include_types)?;
    let detector = parse_detector_json(config.detector_json)?;
    let detector_units = detector_report_units(&detector);
    let gene_by_accession: HashMap<&str, usize> = config
        .index
        .genes
        .iter()
        .enumerate()
        .map(|(idx, gene)| (config.index.string(gene.protein_accession), idx))
        .collect();

    let mut rows = Vec::new();
    for row in &amrfinder_rows {
        rows.push(truth_evidence_row(
            row,
            config.index,
            &gene_by_accession,
            &contigs,
            &detector_units,
            config.min_gene_fraction,
            config.min_family_fraction,
        )?);
    }

    Ok(TruthKmerEvidenceReport {
        assembly: config.assembly_path.display().to_string(),
        detector_json: config.detector_json.display().to_string(),
        amrfinder_tsv: config.amrfinder_tsv.display().to_string(),
        index_k: config.index.k,
        truth_count: rows.len(),
        rows,
    })
}

pub fn debug_amrfinder_misses(config: DebugMissesConfig<'_>) -> anyhow::Result<DebugMissesReport> {
    ensure!(
        config.refinement_k > 0 && config.refinement_k <= config.index.k,
        "refinement k must be between 1 and index k"
    );
    ensure!(
        config.index.alphabet == IndexAlphabet::Dna,
        "AMRFinder miss debugging only supports DNA indexes"
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
        .map(|(idx, gene)| (config.index.string(gene.protein_accession), idx))
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
        let unit_id = gene.report_unit_id as usize;
        let diagnostic = index.unit_specific_kmers(unit_id);
        let matched: HashSet<u64> = diagnostic.intersection(&interval_kmers).copied().collect();
        let missing: Vec<u64> = diagnostic.difference(&interval_kmers).copied().collect();
        let family_diagnostic = hierarchy_unit_kmers(index, index.string(gene.gene_group));
        let family_matched = family_diagnostic.intersection(&interval_kmers).count();
        let lowk_target = refinement_target(&missing, index.k, refinement_k, false);
        let split_target = refinement_target(&missing, index.k, refinement_k, true);
        let interval_lowk = distinct_lowk(&interval, refinement_k);
        let interval_split = distinct_split(&interval, refinement_k);
        (
            Some(index.string(gene.id).to_string()),
            Some(index.string(gene.element_symbol).to_string()),
            Some(index.string(gene.gene_symbol).to_string()),
            Some(index.string(gene.allele_symbol).to_string()),
            Some(index.string(gene.gene_group).to_string()),
            Some(index.string(gene.hierarchy_node).to_string()),
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
    let include_types = HashSet::from(["AMR"]);
    Ok(parse_amrfinder_rows_with_types(path, &include_types)?
        .into_iter()
        .filter(|row| row.subtype == "AMR")
        .collect())
}

fn parse_amrfinder_rows_with_types(
    path: &Path,
    include_types: &HashSet<&str>,
) -> anyhow::Result<Vec<AmrfinderRow>> {
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
        if !include_types.contains(type_value) {
            continue;
        }
        rows.push(AmrfinderRow {
            element_symbol: tsv_field(&fields, &columns, "Element symbol").to_string(),
            hierarchy_node: tsv_field(&fields, &columns, "Hierarchy node").to_string(),
            method: tsv_field(&fields, &columns, "Method").to_string(),
            type_name: type_value.to_string(),
            subtype: tsv_field(&fields, &columns, "Subtype").to_string(),
            class_name: tsv_field(&fields, &columns, "Class").to_string(),
            subclass: tsv_field(&fields, &columns, "Subclass").to_string(),
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

#[derive(Debug, Clone)]
struct UnitEvidence {
    unit_id: usize,
    total: usize,
    matched: usize,
    fraction: f64,
}

fn truth_evidence_row(
    row: &AmrfinderRow,
    index: &AmrIndex,
    gene_by_accession: &HashMap<&str, usize>,
    contigs: &HashMap<String, Vec<u8>>,
    detector_units: &HashSet<String>,
    min_gene_fraction: f64,
    min_family_fraction: f64,
) -> anyhow::Result<TruthKmerEvidenceRow> {
    let interval = match extract_interval(contigs, row) {
        Ok(interval) => interval,
        Err(_) => {
            return Ok(empty_truth_evidence_row(
                row,
                detector_units,
                "interval_unusable",
            ));
        }
    };
    let interval_kmers = distinct_kmers(&interval, index.k);
    let gene_id = gene_by_accession
        .get(row.closest_reference_accession.as_str())
        .copied();
    let mut candidates = Vec::<UnitEvidence>::new();
    let mut exact = UnitEvidence::empty();
    let mut family = UnitEvidence::empty();

    if let Some(gene_id) = gene_id {
        let gene = &index.genes[gene_id];
        exact = unit_evidence(index, gene.report_unit_id as usize, &interval_kmers);
        candidates.push(exact.clone());

        for unit_id in candidate_unit_ids_for_gene(
            index,
            gene.report_unit_id,
            index.string(gene.hierarchy_node),
        ) {
            let evidence = unit_evidence(index, unit_id as usize, &interval_kmers);
            if index.units[unit_id as usize].kind() == ReportUnitKind::HierarchyNode
                && index.string(index.units[unit_id as usize].hierarchy_node)
                    == index.string(gene.hierarchy_node)
            {
                family = evidence.clone();
            }
            candidates.push(evidence);
        }
    } else if !row.hierarchy_node.is_empty() {
        for (unit_id, unit) in index.units.iter().enumerate() {
            if unit.kind() == ReportUnitKind::HierarchyNode
                && index.string(unit.hierarchy_node) == row.hierarchy_node
            {
                let evidence = unit_evidence(index, unit_id, &interval_kmers);
                family = evidence.clone();
                candidates.push(evidence);
            }
        }
    }

    candidates.sort_by(|left, right| {
        right
            .fraction
            .partial_cmp(&left.fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| right.matched.cmp(&left.matched))
    });
    candidates.dedup_by_key(|evidence| evidence.unit_id);
    let best = candidates
        .first()
        .cloned()
        .unwrap_or_else(UnitEvidence::empty);
    let supported = gene_id.is_some() || !candidates.is_empty();
    let covering_detector_units = detector_units_for_truth(row, detector_units, index);
    let covered_by_detector = !covering_detector_units.is_empty();
    let category = recall_failure_category(
        supported,
        interval.is_empty() || interval_kmers.is_empty(),
        covered_by_detector,
        &best,
        min_gene_fraction.min(min_family_fraction),
    );

    Ok(TruthKmerEvidenceRow {
        element_symbol: row.element_symbol.clone(),
        hierarchy_node: row.hierarchy_node.clone(),
        method: row.method.clone(),
        type_name: row.type_name.clone(),
        subtype: row.subtype.clone(),
        class_name: row.class_name.clone(),
        subclass: row.subclass.clone(),
        contig_id: row.contig_id.clone(),
        start: row.start,
        stop: row.stop,
        strand: row.strand.clone(),
        coverage_of_reference: row.coverage_of_reference,
        identity_to_reference: row.identity_to_reference,
        closest_reference_accession: row.closest_reference_accession.clone(),
        covered_by_detector,
        covering_detector_units,
        truth_supported_by_index: supported,
        best_index_unit: best.unit_key(index),
        best_index_unit_type: best.unit_type(index),
        best_index_unit_label: best.unit_label(index),
        best_diagnostic_total: best.total,
        best_diagnostic_matched: best.matched,
        best_diagnostic_missing: best.total.saturating_sub(best.matched),
        best_diagnostic_fraction: best.fraction,
        exact_diagnostic_total: exact.total,
        exact_diagnostic_matched: exact.matched,
        exact_diagnostic_fraction: exact.fraction,
        family_diagnostic_total: family.total,
        family_diagnostic_matched: family.matched,
        family_diagnostic_fraction: family.fraction,
        interval_length: interval.len(),
        interval_distinct_kmers: interval_kmers.len(),
        recall_failure_category: category.to_string(),
    })
}

fn empty_truth_evidence_row(
    row: &AmrfinderRow,
    detector_units: &HashSet<String>,
    category: &str,
) -> TruthKmerEvidenceRow {
    let covering_detector_units = detector_units_for_truth_without_index(row, detector_units);
    TruthKmerEvidenceRow {
        element_symbol: row.element_symbol.clone(),
        hierarchy_node: row.hierarchy_node.clone(),
        method: row.method.clone(),
        type_name: row.type_name.clone(),
        subtype: row.subtype.clone(),
        class_name: row.class_name.clone(),
        subclass: row.subclass.clone(),
        contig_id: row.contig_id.clone(),
        start: row.start,
        stop: row.stop,
        strand: row.strand.clone(),
        coverage_of_reference: row.coverage_of_reference,
        identity_to_reference: row.identity_to_reference,
        closest_reference_accession: row.closest_reference_accession.clone(),
        covered_by_detector: !covering_detector_units.is_empty(),
        covering_detector_units,
        truth_supported_by_index: false,
        best_index_unit: None,
        best_index_unit_type: None,
        best_index_unit_label: None,
        best_diagnostic_total: 0,
        best_diagnostic_matched: 0,
        best_diagnostic_missing: 0,
        best_diagnostic_fraction: 0.0,
        exact_diagnostic_total: 0,
        exact_diagnostic_matched: 0,
        exact_diagnostic_fraction: 0.0,
        family_diagnostic_total: 0,
        family_diagnostic_matched: 0,
        family_diagnostic_fraction: 0.0,
        interval_length: 0,
        interval_distinct_kmers: 0,
        recall_failure_category: category.to_string(),
    }
}

impl UnitEvidence {
    fn empty() -> Self {
        Self {
            unit_id: usize::MAX,
            total: 0,
            matched: 0,
            fraction: 0.0,
        }
    }

    fn unit_key(&self, index: &AmrIndex) -> Option<String> {
        index
            .units
            .get(self.unit_id)
            .map(|unit| index.string(unit.id).to_string())
    }

    fn unit_type(&self, index: &AmrIndex) -> Option<String> {
        index
            .units
            .get(self.unit_id)
            .map(|unit| unit.kind().as_str().to_string())
    }

    fn unit_label(&self, index: &AmrIndex) -> Option<String> {
        index
            .units
            .get(self.unit_id)
            .map(|unit| index.string(unit.label).to_string())
    }
}

fn unit_evidence(index: &AmrIndex, unit_id: usize, interval_kmers: &HashSet<u64>) -> UnitEvidence {
    let diagnostic = index.unit_specific_kmers(unit_id);
    let matched = diagnostic.intersection(interval_kmers).count();
    UnitEvidence {
        unit_id,
        total: diagnostic.len(),
        matched,
        fraction: fraction(matched, diagnostic.len()),
    }
}

fn candidate_unit_ids_for_gene(
    index: &AmrIndex,
    report_unit_id: UnitId,
    hierarchy_node: &str,
) -> Vec<UnitId> {
    let mut ids = Vec::new();
    ids.push(report_unit_id);
    if let Some(unit) = index.units.get(report_unit_id as usize) {
        ids.extend(unit.ancestor_unit_ids.iter().copied());
    }
    ids.extend(index.units.iter().enumerate().filter_map(|(idx, unit)| {
        (unit.kind() == ReportUnitKind::HierarchyNode
            && index.string(unit.hierarchy_node) == hierarchy_node)
            .then_some(idx as UnitId)
    }));
    ids.sort_unstable();
    ids.dedup();
    ids
}

fn detector_report_units(result: &DetectionResult) -> HashSet<String> {
    result.hits.iter().map(detector_report_unit).collect()
}

fn detector_report_unit(hit: &crate::detect::DetectionHit) -> String {
    if hit.call_type == "gene_group" {
        hit.unit_id.clone()
    } else {
        hit.hierarchy_node
            .clone()
            .or_else(|| (!hit.gene_group.is_empty()).then_some(hit.gene_group.clone()))
            .or_else(|| hit.element_symbol.clone())
            .unwrap_or_else(|| hit.unit_id.clone())
    }
}

fn detector_units_for_truth(
    row: &AmrfinderRow,
    detector_units: &HashSet<String>,
    index: &AmrIndex,
) -> Vec<String> {
    let mut units: Vec<String> = detector_units
        .iter()
        .filter(|unit| detector_unit_covers_truth(unit, row, index))
        .cloned()
        .collect();
    units.sort();
    units
}

fn detector_units_for_truth_without_index(
    row: &AmrfinderRow,
    detector_units: &HashSet<String>,
) -> Vec<String> {
    let mut units: Vec<String> = detector_units
        .iter()
        .filter(|unit| *unit == &row.hierarchy_node || *unit == &row.element_symbol)
        .cloned()
        .collect();
    units.sort();
    units
}

fn detector_unit_covers_truth(unit: &str, row: &AmrfinderRow, index: &AmrIndex) -> bool {
    if unit == row.hierarchy_node || unit == row.element_symbol {
        return true;
    }
    let Some(detector_unit) = index
        .units
        .iter()
        .find(|candidate| index.string(candidate.id) == unit)
    else {
        return false;
    };
    if detector_unit.kind() != ReportUnitKind::HierarchyNode {
        return false;
    }
    index.units.iter().any(|truth_unit| {
        index.string(truth_unit.hierarchy_node) == row.hierarchy_node
            && truth_unit
                .ancestor_unit_ids
                .iter()
                .any(|&ancestor_id| index.string(index.units[ancestor_id as usize].id) == unit)
    })
}

fn recall_failure_category(
    supported: bool,
    interval_unusable: bool,
    covered_by_detector: bool,
    best: &UnitEvidence,
    threshold: f64,
) -> &'static str {
    if covered_by_detector {
        return "covered_by_detector";
    }
    if !supported {
        return "unsupported_truth";
    }
    if interval_unusable {
        return "interval_unusable";
    }
    if best.matched == 0 {
        return "no_kmer_evidence";
    }
    if best.fraction >= threshold {
        return "covered_by_kmers_not_reported";
    }
    if best.fraction >= threshold * 0.5 {
        return "near_threshold";
    }
    "weak_kmer_evidence"
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

fn hierarchy_unit_kmers(index: &AmrIndex, hierarchy_node: &str) -> HashSet<u64> {
    let Some(unit_id) = index.units.iter().position(|unit| {
        unit.kind() == ReportUnitKind::HierarchyNode
            && index.string(unit.hierarchy_node) == hierarchy_node
    }) else {
        return HashSet::new();
    };
    index.unit_specific_kmers(unit_id)
}

fn reference_kmers_by_accession(
    db_dir: &Path,
    k: usize,
) -> anyhow::Result<HashMap<String, HashSet<u64>>> {
    let mut by_accession = HashMap::<String, HashSet<u64>>::new();
    for reference in load_amrfinder_references(db_dir, &[ReferenceType::Amr])? {
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
