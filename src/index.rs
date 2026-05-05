use crate::amrfinder_db::AmrReference;
use crate::kmer::DnaKmerIter;
use anyhow::{Context, ensure};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

const INDEX_MAGIC: &[u8; 8] = b"SHAMR001";
const ASSIGN_GENE: u8 = 1;
const ASSIGN_FAMILY: u8 = 2;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneEntry {
    pub id: String,
    pub element_symbol: String,
    pub gene_symbol: String,
    pub allele_symbol: String,
    pub protein_accession: String,
    pub nucleotide_accession: String,
    pub family: String,
    pub class_name: String,
    pub subclass: String,
    pub hierarchy_node: String,
    pub product: String,
    pub length: usize,
    pub gene_specific_kmers: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmrIndex {
    pub magic: [u8; 8],
    pub db_version: String,
    pub k: usize,
    pub genes: Vec<GeneEntry>,
    pub families: Vec<String>,
    pub family_specific_kmers: Vec<usize>,
    pub kmer_codes: Vec<u64>,
    pub assignment_kind: Vec<u8>,
    pub assignment_id: Vec<u32>,
    pub gene_kmer_offsets: Vec<u32>,
    pub gene_kmers: Vec<u64>,
}

#[derive(Debug, Clone)]
pub struct IndexBuildConfig {
    pub k: usize,
}

impl Default for IndexBuildConfig {
    fn default() -> Self {
        Self { k: 31 }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KmerAssignment {
    Gene(usize),
    Family(usize),
}

impl AmrIndex {
    pub fn lookup(&self, code: u64) -> Option<KmerAssignment> {
        let idx = self.kmer_codes.binary_search(&code).ok()?;
        match self.assignment_kind[idx] {
            ASSIGN_GENE => Some(KmerAssignment::Gene(self.assignment_id[idx] as usize)),
            ASSIGN_FAMILY => Some(KmerAssignment::Family(self.assignment_id[idx] as usize)),
            _ => None,
        }
    }

    pub fn gene_specific_kmers(&self, gene_id: usize) -> &[u64] {
        let start = self.gene_kmer_offsets[gene_id] as usize;
        let end = self.gene_kmer_offsets[gene_id + 1] as usize;
        &self.gene_kmers[start..end]
    }

    pub fn stats_string(&self) -> String {
        let genes_without_gene_specific_kmers = self
            .genes
            .iter()
            .filter(|gene| gene.gene_specific_kmers == 0)
            .count();
        format!(
            "db_version={}\nk={}\ngenes={}\ngenes_without_gene_specific_kmers={}\nfamilies={}\nretained_kmers={}\ngene_specific_kmers={}\nfamily_specific_kmers={}\nindex_bytes_uncompressed_estimate={}\n",
            self.db_version,
            self.k,
            self.genes.len(),
            genes_without_gene_specific_kmers,
            self.families.len(),
            self.kmer_codes.len(),
            self.gene_kmers.len(),
            self.family_specific_kmers.iter().sum::<usize>(),
            self.estimated_bytes(),
        )
    }

    fn estimated_bytes(&self) -> usize {
        self.kmer_codes.len() * 8
            + self.assignment_kind.len()
            + self.assignment_id.len() * 4
            + self.gene_kmers.len() * 8
            + self.gene_kmer_offsets.len() * 4
    }
}

pub fn build_index(
    references: &[AmrReference],
    config: &IndexBuildConfig,
) -> anyhow::Result<AmrIndex> {
    ensure!(
        (1..=31).contains(&config.k),
        "DNA k must be between 1 and 31"
    );

    let mut family_to_id = HashMap::<String, usize>::new();
    let mut families = Vec::<String>::new();
    let mut genes = Vec::<GeneEntry>::new();
    let mut raw_kmers = HashMap::<u64, Vec<(usize, usize)>>::new();

    for reference in references {
        let family_id = *family_to_id
            .entry(reference.family.clone())
            .or_insert_with(|| {
                families.push(reference.family.clone());
                families.len() - 1
            });
        let gene_id = genes.len();
        let gene_label = if reference.allele_symbol.is_empty() {
            reference.gene_symbol.clone()
        } else {
            reference.allele_symbol.clone()
        };
        genes.push(GeneEntry {
            id: format!("{}|{}", gene_label, reference.protein_accession),
            element_symbol: reference.element_symbol.clone(),
            gene_symbol: reference.gene_symbol.clone(),
            allele_symbol: reference.allele_symbol.clone(),
            protein_accession: reference.protein_accession.clone(),
            nucleotide_accession: reference.nucleotide_accession.clone(),
            family: reference.family.clone(),
            class_name: reference.class_name.clone(),
            subclass: reference.subclass.clone(),
            hierarchy_node: reference.hierarchy_node.clone(),
            product: reference.product.clone(),
            length: reference.seq.len(),
            gene_specific_kmers: 0,
        });

        let Some(iter) = DnaKmerIter::new(&reference.seq, config.k) else {
            continue;
        };
        let mut seen_in_gene = HashSet::new();
        for (_pos, kmer) in iter {
            if seen_in_gene.insert(kmer) {
                raw_kmers
                    .entry(kmer)
                    .or_default()
                    .push((gene_id, family_id));
            }
        }
    }

    let mut retained = Vec::<(u64, u8, u32)>::new();
    let mut per_gene_kmers = vec![Vec::<u64>::new(); genes.len()];
    let mut family_specific_kmers = vec![0usize; families.len()];

    for (kmer, refs) in raw_kmers {
        let gene_ids: HashSet<usize> = refs.iter().map(|(gene_id, _)| *gene_id).collect();
        if gene_ids.len() == 1 {
            let gene_id = *gene_ids.iter().next().unwrap();
            genes[gene_id].gene_specific_kmers += 1;
            per_gene_kmers[gene_id].push(kmer);
            retained.push((kmer, ASSIGN_GENE, gene_id as u32));
            continue;
        }

        let family_ids: HashSet<usize> = refs.iter().map(|(_, family_id)| *family_id).collect();
        if family_ids.len() == 1 {
            let family_id = *family_ids.iter().next().unwrap();
            family_specific_kmers[family_id] += 1;
            retained.push((kmer, ASSIGN_FAMILY, family_id as u32));
        }
    }

    retained.sort_by_key(|(kmer, _, _)| *kmer);
    for kmers in &mut per_gene_kmers {
        kmers.sort_unstable();
    }

    let mut gene_kmer_offsets = Vec::with_capacity(genes.len() + 1);
    let mut gene_kmers = Vec::new();
    gene_kmer_offsets.push(0);
    for kmers in per_gene_kmers {
        gene_kmers.extend(kmers);
        gene_kmer_offsets.push(gene_kmers.len() as u32);
    }

    let db_version = references
        .first()
        .map(|reference| reference.db_version.clone())
        .unwrap_or_else(|| "unknown".to_string());

    Ok(AmrIndex {
        magic: *INDEX_MAGIC,
        db_version,
        k: config.k,
        genes,
        families,
        family_specific_kmers,
        kmer_codes: retained.iter().map(|(kmer, _, _)| *kmer).collect(),
        assignment_kind: retained.iter().map(|(_, kind, _)| *kind).collect(),
        assignment_id: retained.iter().map(|(_, _, id)| *id).collect(),
        gene_kmer_offsets,
        gene_kmers,
    })
}

pub fn save_index(index: &AmrIndex, path: &Path) -> anyhow::Result<()> {
    let bytes = bincode::serialize(index).context("serialize AMR index")?;
    fs::write(path, bytes).with_context(|| format!("write index {}", path.display()))
}

pub fn load_index(path: &Path) -> anyhow::Result<AmrIndex> {
    let bytes = fs::read(path).with_context(|| format!("read index {}", path.display()))?;
    let index: AmrIndex = bincode::deserialize(&bytes).context("deserialize AMR index")?;
    ensure!(index.magic == *INDEX_MAGIC, "unsupported AMR index format");
    ensure!(
        index.kmer_codes.len() == index.assignment_kind.len()
            && index.kmer_codes.len() == index.assignment_id.len(),
        "corrupt AMR index assignment arrays"
    );
    ensure!(
        index.gene_kmer_offsets.len() == index.genes.len() + 1,
        "corrupt AMR index gene offsets"
    );
    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference(id: &str, family: &str, seq: &[u8]) -> AmrReference {
        AmrReference {
            protein_accession: id.to_string(),
            nucleotide_accession: String::new(),
            element_symbol: id.to_string(),
            gene_symbol: id.to_string(),
            allele_symbol: id.to_string(),
            product: String::new(),
            family: family.to_string(),
            class_name: String::new(),
            subclass: String::new(),
            hierarchy_node: String::new(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            db_version: "test".to_string(),
            seq: seq.to_vec(),
        }
    }

    #[test]
    fn builds_lookup_arrays() {
        let refs = vec![
            reference("g1", "f1", b"ACGTACGTAC"),
            reference("g2", "f2", b"TTTTTCCCCCA"),
        ];
        let index = build_index(&refs, &IndexBuildConfig { k: 5 }).unwrap();
        assert_eq!(index.genes.len(), 2);
        assert!(!index.kmer_codes.is_empty());
        assert_eq!(index.gene_kmer_offsets.len(), 3);
    }
}
