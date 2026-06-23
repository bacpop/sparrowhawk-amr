use crate::fasta::parse_fasta_bytes;
use crate::index::{AmrIndex, IndexAlphabet, ReportUnit, ReportUnitKind};
use crate::kmer::{DnaKmerIter, ProteinKmerIter};
use anyhow::ensure;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum QueryKind {
    Direct,
    Cds,
    ProteinCds,
}

/// TEST This temporary enum was done for doing tests and seeing if getting lower k-vals
/// after first matches, or using split k-mers might help refine the results
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RefinementMode {
    None,
    Split,
    LowK,
}

impl std::fmt::Display for RefinementMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, "none"),
            Self::Split => write!(f, "split"),
            Self::LowK => write!(f, "lowk"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectParams {
    pub min_gene_fraction: f64,
    #[serde(alias = "min_gene_group_fraction")]
    pub min_gene_group_fraction: f64,
    pub seed_gene_fraction: f64,
    pub seed_gene_hits: usize,
    pub refinement_mode: RefinementMode,
    pub refinement_k: usize,
}

impl Default for DetectParams {
    fn default() -> Self {
        Self {
            min_gene_fraction: 0.10,
            min_gene_group_fraction: 0.10,
            seed_gene_fraction: 0.01,
            seed_gene_hits: 3,
            refinement_mode: RefinementMode::None,
            refinement_k: 21,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectionHit {
    pub query_id: String,
    pub query_kind: QueryKind,
    pub unit_id: String,
    pub unit_label: String,
    pub gene_id: Option<String>,
    pub element_symbol: Option<String>,
    pub gene_symbol: Option<String>,
    pub allele_symbol: Option<String>,
    #[serde(alias = "family")]
    pub gene_group: String,
    pub hierarchy_node: Option<String>,
    pub class_name: Option<String>,
    pub subclass: Option<String>,
    pub type_name: Option<String>,
    pub subtype: Option<String>,
    pub member_count: usize,
    pub start: usize,
    pub end: usize,
    pub call_stage: String,
    pub first_pass_distinct: usize,
    pub first_pass_total: usize,
    pub first_pass_diagnostic_total: usize,
    pub first_pass_fraction: f64,
    pub refinement_distinct: usize,
    pub refinement_total: usize,
    pub refinement_diagnostic_total: usize,
    pub refinement_fraction: f64,
    pub call_fraction: f64,
    pub call_type: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectionResult {
    pub sample_name: String,
    pub database_version: String,
    pub query_kind: QueryKind,
    #[serde(default)]
    pub index_alphabet: IndexAlphabet,
    pub index_k: usize,
    pub refinement_mode: RefinementMode,
    pub refinement_k: usize,
    pub hits: Vec<DetectionHit>,
    pub gene_count: usize,
    pub gene_group_count: usize,
}

#[derive(Debug, Clone, Default)]
struct HitAccumulator {
    count: usize,
    distinct: HashSet<u64>,
    min_pos: usize,
    max_pos: usize,
}

impl HitAccumulator {
    fn add(&mut self, kmer: u64, pos: usize, k: usize) {
        self.count += 1;
        self.distinct.insert(kmer);
        if self.count == 1 {
            self.min_pos = pos;
            self.max_pos = pos + k;
        } else {
            self.min_pos = self.min_pos.min(pos);
            self.max_pos = self.max_pos.max(pos + k);
        }
    }
}

pub fn detect_fasta(
    index: &AmrIndex,
    fasta_bytes: &[u8],
    sample_name: &str,
    query_kind: QueryKind,
    params: &DetectParams,
) -> anyhow::Result<DetectionResult> {
    ensure!(
        index.alphabet == expected_alphabet(query_kind),
        "{} index cannot be used for {:?} detection",
        index.alphabet.as_str(),
        query_kind
    );
    let records = parse_fasta_bytes(fasta_bytes)?;
    let mut hits = Vec::new();

    for record in records {
        let mut unit_hits = HashMap::<usize, HitAccumulator>::new();

        match index.alphabet {
            IndexAlphabet::Dna => {
                if let Some(iter) = DnaKmerIter::new(&record.seq, index.k) {
                    for (pos, kmer) in iter {
                        if let Some(unit_id) = index.lookup(kmer) {
                            unit_hits
                                .entry(unit_id)
                                .or_default()
                                .add(kmer, pos, index.k);
                        }
                    }
                }
            }
            IndexAlphabet::Protein => {
                if let Some(iter) = ProteinKmerIter::new(&record.seq, index.k) {
                    for (pos, kmer) in iter {
                        if let Some(unit_id) = index.lookup(kmer) {
                            unit_hits
                                .entry(unit_id)
                                .or_default()
                                .add(kmer, pos, index.k);
                        }
                    }
                }
            }
        }

        let mut suppressed_hierarchy_units = HashSet::<usize>::new();
        for (&unit_id, acc) in &unit_hits {
            let unit = &index.units[unit_id];
            if unit.kind() != ReportUnitKind::ExactGene || unit.diagnostic_kmers == 0 {
                continue;
            }
            let fraction = acc.distinct.len() as f64 / unit.diagnostic_kmers as f64;
            if fraction < params.min_gene_fraction {
                continue;
            }
            suppressed_hierarchy_units.extend(
                unit.ancestor_unit_ids
                    .iter()
                    .map(|&ancestor_id| ancestor_id as usize),
            );
            hits.push(unit_hit(
                index, unit, &record.id, query_kind, index.k, acc, fraction,
            ));
        }

        let mut hierarchy_candidates: Vec<(usize, &HitAccumulator, f64)> = unit_hits
            .iter()
            .filter_map(|(&unit_id, acc)| {
                let unit = &index.units[unit_id];
                if unit.kind() != ReportUnitKind::HierarchyNode || unit.diagnostic_kmers == 0 {
                    return None;
                }
                let fraction = acc.distinct.len() as f64 / unit.diagnostic_kmers as f64;
                (fraction >= params.min_gene_group_fraction).then_some((unit_id, acc, fraction))
            })
            .collect();
        hierarchy_candidates.sort_by_key(|(unit_id, _, _)| index.units[*unit_id].member_count);

        for (unit_id, acc, fraction) in hierarchy_candidates {
            if suppressed_hierarchy_units.contains(&unit_id) {
                continue;
            }
            let unit = &index.units[unit_id];
            suppressed_hierarchy_units.extend(
                unit.ancestor_unit_ids
                    .iter()
                    .map(|&ancestor_id| ancestor_id as usize),
            );
            hits.push(unit_hit(
                index, unit, &record.id, query_kind, index.k, acc, fraction,
            ));
        }
    }

    let gene_count = hits.iter().filter(|hit| hit.call_type == "gene").count();
    let gene_group_count = hits
        .iter()
        .filter(|hit| hit.call_type == "gene_group")
        .count();
    Ok(DetectionResult {
        sample_name: sample_name.to_string(),
        database_version: index.db_version.clone(),
        query_kind,
        index_alphabet: index.alphabet,
        index_k: index.k,
        refinement_mode: params.refinement_mode,
        refinement_k: params.refinement_k,
        hits,
        gene_count,
        gene_group_count,
    })
}

pub fn detect_protein_fasta(
    index: &AmrIndex,
    fasta_bytes: &[u8],
    sample_name: &str,
    params: &DetectParams,
) -> anyhow::Result<DetectionResult> {
    detect_fasta(
        index,
        fasta_bytes,
        sample_name,
        QueryKind::ProteinCds,
        params,
    )
}

fn expected_alphabet(query_kind: QueryKind) -> IndexAlphabet {
    match query_kind {
        QueryKind::Direct | QueryKind::Cds => IndexAlphabet::Dna,
        QueryKind::ProteinCds => IndexAlphabet::Protein,
    }
}


// Recover all the info, including metadata, from the index
fn unit_hit(
    index: &AmrIndex,
    unit: &ReportUnit,
    query_id: &str,
    query_kind: QueryKind,
    k: usize,
    acc: &HitAccumulator,
    call_fraction: f64,
) -> DetectionHit {
    let first_fraction = if unit.diagnostic_kmers == 0 {
        0.0
    } else {
        acc.distinct.len() as f64 / unit.diagnostic_kmers as f64
    };

    DetectionHit {
        query_id: query_id.to_string(),
        query_kind,
        unit_id: index.string(unit.id).to_string(),
        unit_label: index.string(unit.label).to_string(),
        gene_id: unit.gene_id.map(|_| index.string(unit.id).to_string()),
        element_symbol: index.optional_string(unit.element_symbol),
        gene_symbol: index.optional_string(unit.gene_symbol),
        allele_symbol: index.optional_string(unit.allele_symbol),
        gene_group: index.string(unit.gene_group).to_string(),
        hierarchy_node: non_empty_option(index.string(unit.hierarchy_node)),
        class_name: non_empty_option(index.string(unit.class_name)),
        subclass: non_empty_option(index.string(unit.subclass)),
        type_name: non_empty_option(index.string(unit.type_name)),
        subtype: non_empty_option(index.string(unit.subtype)),
        member_count: unit.member_count,
        start: acc.min_pos,
        end: acc.max_pos,
        call_stage: format!("k{k}"),
        first_pass_distinct: acc.distinct.len(),
        first_pass_total: acc.count,
        first_pass_diagnostic_total: unit.diagnostic_kmers,
        first_pass_fraction: first_fraction,
        refinement_distinct: 0,
        refinement_total: 0,
        refinement_diagnostic_total: 0,
        refinement_fraction: 0.0,
        call_fraction,
        call_type: unit.call_type().to_string(),
    }
}

fn non_empty_option(value: &str) -> Option<String> {
    (!value.is_empty()).then(|| value.to_string())
}





// =============================== TESTS

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amrfinder_db::{AmrReference, HierarchyNode};
    use crate::index::{IndexBuildConfig, build_index};

    #[test]
    fn detects_simple_gene() {
        let refs = vec![AmrReference {
            protein_accession: "p1".to_string(),
            nucleotide_accession: "n1".to_string(),
            element_symbol: "geneA".to_string(),
            gene_symbol: "geneA".to_string(),
            allele_symbol: "geneA".to_string(),
            product: String::new(),
            family: "famA".to_string(),
            class_name: "CLASS".to_string(),
            subclass: "SUB".to_string(),
            hierarchy_node: "node".to_string(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            hierarchy_path: vec![HierarchyNode {
                node_id: "node".to_string(),
                parent_node_id: String::new(),
                symbol: "node".to_string(),
                class_name: "CLASS".to_string(),
                subclass: "SUB".to_string(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: 2,
            }],
            db_version: "test".to_string(),
            seq: b"ACGTACGTACGT".to_vec(),
        }];
        let index = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Dna,
                k: 5,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        let params = DetectParams {
            min_gene_fraction: 0.5,
            ..DetectParams::default()
        };
        let result = detect_fasta(
            &index,
            b">contig\nTTTACGTACGTACGTTTT\n",
            "sample",
            QueryKind::Direct,
            &params,
        )
        .unwrap();
        assert_eq!(result.gene_count, 1);
        assert_eq!(result.hits[0].type_name.as_deref(), Some("AMR"));
        assert_eq!(result.hits[0].subtype.as_deref(), Some("AMR"));
    }

    #[test]
    fn detects_cds_query_kind_against_dna_index() {
        let refs = vec![AmrReference {
            protein_accession: "p1".to_string(),
            nucleotide_accession: "n1".to_string(),
            element_symbol: "geneA".to_string(),
            gene_symbol: "geneA".to_string(),
            allele_symbol: "geneA".to_string(),
            product: String::new(),
            family: "famA".to_string(),
            class_name: "CLASS".to_string(),
            subclass: "SUB".to_string(),
            hierarchy_node: "node".to_string(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            hierarchy_path: vec![HierarchyNode {
                node_id: "node".to_string(),
                parent_node_id: String::new(),
                symbol: "node".to_string(),
                class_name: "CLASS".to_string(),
                subclass: "SUB".to_string(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: 2,
            }],
            db_version: "test".to_string(),
            seq: b"ACGTACGTACGT".to_vec(),
        }];
        let index = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Dna,
                k: 5,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        let result = detect_fasta(
            &index,
            b">gene_1
ACGTACGTACGT
",
            "sample",
            QueryKind::Cds,
            &DetectParams {
                min_gene_fraction: 0.5,
                ..DetectParams::default()
            },
        )
        .unwrap();
        assert_eq!(result.query_kind, QueryKind::Cds);
        assert_eq!(result.hits[0].query_id, "gene_1");
    }

    #[test]
    fn rejects_protein_query_against_dna_index() {
        let refs = vec![AmrReference {
            protein_accession: "p1".to_string(),
            nucleotide_accession: "n1".to_string(),
            element_symbol: "geneA".to_string(),
            gene_symbol: "geneA".to_string(),
            allele_symbol: "geneA".to_string(),
            product: String::new(),
            family: "famA".to_string(),
            class_name: "CLASS".to_string(),
            subclass: "SUB".to_string(),
            hierarchy_node: "node".to_string(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            hierarchy_path: vec![HierarchyNode {
                node_id: "node".to_string(),
                parent_node_id: String::new(),
                symbol: "node".to_string(),
                class_name: "CLASS".to_string(),
                subclass: "SUB".to_string(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: 2,
            }],
            db_version: "test".to_string(),
            seq: b"ACGTACGTACGT".to_vec(),
        }];
        let index = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Dna,
                k: 5,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        let err = detect_protein_fasta(
            &index,
            b">protein\nMKTAA\n",
            "sample",
            &DetectParams::default(),
        )
        .unwrap_err();
        assert!(err.to_string().contains("dna index"));
    }

    #[test]
    fn detects_simple_protein_gene() {
        let refs = vec![AmrReference {
            protein_accession: "p1".to_string(),
            nucleotide_accession: "n1".to_string(),
            element_symbol: "geneA".to_string(),
            gene_symbol: "geneA".to_string(),
            allele_symbol: "geneA".to_string(),
            product: String::new(),
            family: "famA".to_string(),
            class_name: "CLASS".to_string(),
            subclass: "SUB".to_string(),
            hierarchy_node: "node".to_string(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            hierarchy_path: vec![HierarchyNode {
                node_id: "node".to_string(),
                parent_node_id: String::new(),
                symbol: "node".to_string(),
                class_name: "CLASS".to_string(),
                subclass: "SUB".to_string(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: 2,
            }],
            db_version: "test".to_string(),
            seq: b"MKTAA".to_vec(),
        }];
        let index = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Protein,
                k: 3,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        let result = detect_protein_fasta(
            &index,
            b">protein\nXXMKTAA\n",
            "sample",
            &DetectParams {
                min_gene_fraction: 0.5,
                ..DetectParams::default()
            },
        )
        .unwrap();
        assert_eq!(result.index_alphabet, IndexAlphabet::Protein);
        assert_eq!(result.gene_count, 1);
    }

    #[test]
    fn rejects_direct_query_against_protein_index() {
        let refs = vec![AmrReference {
            protein_accession: "p1".to_string(),
            nucleotide_accession: "n1".to_string(),
            element_symbol: "geneA".to_string(),
            gene_symbol: "geneA".to_string(),
            allele_symbol: "geneA".to_string(),
            product: String::new(),
            family: "famA".to_string(),
            class_name: "CLASS".to_string(),
            subclass: "SUB".to_string(),
            hierarchy_node: "node".to_string(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            hierarchy_path: vec![HierarchyNode {
                node_id: "node".to_string(),
                parent_node_id: String::new(),
                symbol: "node".to_string(),
                class_name: "CLASS".to_string(),
                subclass: "SUB".to_string(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: 2,
            }],
            db_version: "test".to_string(),
            seq: b"MKTAA".to_vec(),
        }];
        let index = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Protein,
                k: 3,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        let err = detect_fasta(
            &index,
            b">contig\nATGAAAACCGCC\n",
            "sample",
            QueryKind::Direct,
            &DetectParams::default(),
        )
        .unwrap_err();
        assert!(err.to_string().contains("protein index"));
    }
}
