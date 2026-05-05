use crate::fasta::parse_fasta_bytes;
use crate::index::{AmrIndex, KmerAssignment};
use crate::kmer::{DnaKmerIter, SplitKmerIter, decode_kmer};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum QueryKind {
    Direct,
    Cds,
}

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
    pub min_family_fraction: f64,
    pub seed_gene_fraction: f64,
    pub seed_gene_hits: usize,
    pub refinement_mode: RefinementMode,
    pub refinement_k: usize,
}

impl Default for DetectParams {
    fn default() -> Self {
        Self {
            min_gene_fraction: 0.10,
            min_family_fraction: 0.10,
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
    pub gene_id: Option<String>,
    pub element_symbol: Option<String>,
    pub gene_symbol: Option<String>,
    pub allele_symbol: Option<String>,
    pub family: String,
    pub class_name: Option<String>,
    pub subclass: Option<String>,
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
    pub index_k: usize,
    pub refinement_mode: RefinementMode,
    pub refinement_k: usize,
    pub hits: Vec<DetectionHit>,
    pub gene_count: usize,
    pub family_count: usize,
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
    let records = parse_fasta_bytes(fasta_bytes)?;
    let mut hits = Vec::new();

    for record in records {
        let mut gene_hits = HashMap::<usize, HitAccumulator>::new();
        let mut family_hits = HashMap::<usize, HitAccumulator>::new();

        if let Some(iter) = DnaKmerIter::new(&record.seq, index.k) {
            for (pos, kmer) in iter {
                match index.lookup(kmer) {
                    Some(KmerAssignment::Gene(gene_id)) => {
                        gene_hits
                            .entry(gene_id)
                            .or_default()
                            .add(kmer, pos, index.k);
                    }
                    Some(KmerAssignment::Family(family_id)) => {
                        family_hits
                            .entry(family_id)
                            .or_default()
                            .add(kmer, pos, index.k);
                    }
                    None => {}
                }
            }
        }

        let mut called_genes = HashSet::<usize>::new();
        let mut claimed_families = HashSet::<String>::new();
        for (&gene_id, acc) in &gene_hits {
            let gene = &index.genes[gene_id];
            if gene.gene_specific_kmers == 0 {
                continue;
            }
            let fraction = acc.distinct.len() as f64 / gene.gene_specific_kmers as f64;
            if fraction < params.min_gene_fraction {
                continue;
            }
            called_genes.insert(gene_id);
            claimed_families.insert(gene.family.clone());
            hits.push(gene_hit(
                index, gene_id, &record.id, query_kind, acc, "k31", 0, 0, 0, 0.0, fraction,
            ));
        }

        if params.refinement_mode != RefinementMode::None {
            for (&gene_id, acc) in &gene_hits {
                if called_genes.contains(&gene_id) {
                    continue;
                }
                let gene = &index.genes[gene_id];
                if gene.gene_specific_kmers == 0 {
                    continue;
                }
                let first_fraction = acc.distinct.len() as f64 / gene.gene_specific_kmers as f64;
                if acc.distinct.len() < params.seed_gene_hits
                    && first_fraction < params.seed_gene_fraction
                {
                    continue;
                }
                let Some(refined) = refine_gene(index, gene_id, &record.seq, &acc.distinct, params)
                else {
                    continue;
                };
                if refined.diagnostic_total == 0 {
                    continue;
                }
                let refinement_fraction =
                    refined.distinct.len() as f64 / refined.diagnostic_total as f64;
                if refinement_fraction < params.min_gene_fraction {
                    continue;
                }
                called_genes.insert(gene_id);
                claimed_families.insert(gene.family.clone());
                hits.push(gene_hit(
                    index,
                    gene_id,
                    &record.id,
                    query_kind,
                    acc,
                    &params.refinement_mode.to_string(),
                    refined.distinct.len(),
                    refined.count,
                    refined.diagnostic_total,
                    refinement_fraction,
                    refinement_fraction.max(first_fraction),
                ));
            }
        }

        for (&family_id, acc) in &family_hits {
            let family = &index.families[family_id];
            let diagnostic_total = index.family_specific_kmers[family_id];
            if diagnostic_total == 0 || claimed_families.contains(family) {
                continue;
            }
            let fraction = acc.distinct.len() as f64 / diagnostic_total as f64;
            if fraction < params.min_family_fraction {
                continue;
            }
            hits.push(DetectionHit {
                query_id: record.id.clone(),
                query_kind,
                gene_id: None,
                element_symbol: None,
                gene_symbol: None,
                allele_symbol: None,
                family: family.clone(),
                class_name: None,
                subclass: None,
                start: acc.min_pos,
                end: acc.max_pos,
                call_stage: "k31".to_string(),
                first_pass_distinct: acc.distinct.len(),
                first_pass_total: acc.count,
                first_pass_diagnostic_total: diagnostic_total,
                first_pass_fraction: fraction,
                refinement_distinct: 0,
                refinement_total: 0,
                refinement_diagnostic_total: 0,
                refinement_fraction: 0.0,
                call_fraction: fraction,
                call_type: "family".to_string(),
            });
        }
    }

    let gene_count = hits.iter().filter(|hit| hit.call_type == "gene").count();
    let family_count = hits.iter().filter(|hit| hit.call_type == "family").count();
    Ok(DetectionResult {
        sample_name: sample_name.to_string(),
        database_version: index.db_version.clone(),
        query_kind,
        index_k: index.k,
        refinement_mode: params.refinement_mode,
        refinement_k: params.refinement_k,
        hits,
        gene_count,
        family_count,
    })
}

fn gene_hit(
    index: &AmrIndex,
    gene_id: usize,
    query_id: &str,
    query_kind: QueryKind,
    acc: &HitAccumulator,
    stage: &str,
    refinement_distinct: usize,
    refinement_total: usize,
    refinement_diagnostic_total: usize,
    refinement_fraction: f64,
    call_fraction: f64,
) -> DetectionHit {
    let gene = &index.genes[gene_id];
    let first_fraction = if gene.gene_specific_kmers == 0 {
        0.0
    } else {
        acc.distinct.len() as f64 / gene.gene_specific_kmers as f64
    };
    DetectionHit {
        query_id: query_id.to_string(),
        query_kind,
        gene_id: Some(gene.id.clone()),
        element_symbol: Some(gene.element_symbol.clone()),
        gene_symbol: Some(gene.gene_symbol.clone()),
        allele_symbol: Some(gene.allele_symbol.clone()),
        family: gene.family.clone(),
        class_name: Some(gene.class_name.clone()),
        subclass: Some(gene.subclass.clone()),
        start: acc.min_pos,
        end: acc.max_pos,
        call_stage: stage.to_string(),
        first_pass_distinct: acc.distinct.len(),
        first_pass_total: acc.count,
        first_pass_diagnostic_total: gene.gene_specific_kmers,
        first_pass_fraction: first_fraction,
        refinement_distinct,
        refinement_total,
        refinement_diagnostic_total,
        refinement_fraction,
        call_fraction,
        call_type: "gene".to_string(),
    }
}

#[derive(Debug)]
struct RefinementAccumulator {
    count: usize,
    distinct: HashSet<u64>,
    diagnostic_total: usize,
}

fn refine_gene(
    index: &AmrIndex,
    gene_id: usize,
    query_seq: &[u8],
    matched_first_pass: &HashSet<u64>,
    params: &DetectParams,
) -> Option<RefinementAccumulator> {
    let mut target = HashSet::<u64>::new();
    for &kmer in index.gene_specific_kmers(gene_id) {
        if matched_first_pass.contains(&kmer) {
            continue;
        }
        let decoded = decode_kmer(kmer, index.k);
        match params.refinement_mode {
            RefinementMode::None => {}
            RefinementMode::Split => {
                if let Some(iter) = SplitKmerIter::new(&decoded, params.refinement_k) {
                    target.extend(iter.map(|(_, code)| code));
                }
            }
            RefinementMode::LowK => {
                if let Some(iter) = DnaKmerIter::new(&decoded, params.refinement_k) {
                    target.extend(iter.map(|(_, code)| code));
                }
            }
        }
    }
    if target.is_empty() {
        return None;
    }

    let mut distinct = HashSet::new();
    let mut count = 0usize;
    match params.refinement_mode {
        RefinementMode::None => {}
        RefinementMode::Split => {
            if let Some(iter) = SplitKmerIter::new(query_seq, params.refinement_k) {
                for (_pos, code) in iter {
                    if target.contains(&code) {
                        count += 1;
                        distinct.insert(code);
                    }
                }
            }
        }
        RefinementMode::LowK => {
            if let Some(iter) = DnaKmerIter::new(query_seq, params.refinement_k) {
                for (_pos, code) in iter {
                    if target.contains(&code) {
                        count += 1;
                        distinct.insert(code);
                    }
                }
            }
        }
    }

    Some(RefinementAccumulator {
        count,
        distinct,
        diagnostic_total: target.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amrfinder_db::AmrReference;
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
            db_version: "test".to_string(),
            seq: b"ACGTACGTACGT".to_vec(),
        }];
        let index = build_index(&refs, &IndexBuildConfig { k: 5 }).unwrap();
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
    }
}
