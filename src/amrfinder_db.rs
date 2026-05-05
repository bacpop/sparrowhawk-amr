use crate::fasta::{FastaRecord, read_fasta};
use anyhow::Context;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct AmrReference {
    pub protein_accession: String,
    pub nucleotide_accession: String,
    pub element_symbol: String,
    pub gene_symbol: String,
    pub allele_symbol: String,
    pub product: String,
    pub family: String,
    pub class_name: String,
    pub subclass: String,
    pub hierarchy_node: String,
    pub scope: String,
    pub type_name: String,
    pub subtype: String,
    pub reportable: u8,
    pub db_version: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone, Default)]
struct Metadata {
    element_symbol: String,
    class_name: String,
    subclass: String,
    hierarchy_node: String,
    scope: String,
    type_name: String,
    subtype: String,
    reportable: u8,
}

#[derive(Debug, Clone)]
struct CdsHeader {
    protein_accession: String,
    nucleotide_accession: String,
    gene_symbol: String,
    allele_symbol: String,
    product: String,
}

pub fn load_amrfinder_references(db_dir: &Path) -> anyhow::Result<Vec<AmrReference>> {
    let version = fs::read_to_string(db_dir.join("version.txt"))
        .unwrap_or_else(|_| "unknown".to_string())
        .trim()
        .to_string();
    let metadata = load_metadata(db_dir)?;
    let fasta_path = db_dir.join("AMR_CDS.fa");
    let records = read_fasta(&fasta_path)
        .with_context(|| format!("read AMRFinderPlus CDS FASTA {}", fasta_path.display()))?;

    let mut references = Vec::new();
    for record in records {
        let header = parse_cds_header(&record);
        let meta = metadata
            .get(&header.protein_accession)
            .or_else(|| metadata.get(&header.gene_symbol))
            .or_else(|| metadata.get(&header.allele_symbol))
            .cloned()
            .unwrap_or_default();

        if !meta.type_name.is_empty() && meta.type_name != "AMR" {
            continue;
        }
        if !meta.subtype.is_empty() && meta.subtype != "AMR" {
            continue;
        }

        references.push(AmrReference {
            protein_accession: header.protein_accession,
            nucleotide_accession: header.nucleotide_accession,
            element_symbol: fallback_element_symbol(
                &meta.element_symbol,
                &header.gene_symbol,
                &header.allele_symbol,
            ),
            gene_symbol: header.gene_symbol.clone(),
            allele_symbol: header.allele_symbol.clone(),
            product: header.product,
            family: fallback_family(
                &meta.hierarchy_node,
                &header.gene_symbol,
                &header.allele_symbol,
            ),
            class_name: meta.class_name,
            subclass: meta.subclass,
            hierarchy_node: meta.hierarchy_node,
            scope: meta.scope,
            type_name: if meta.type_name.is_empty() {
                "AMR".to_string()
            } else {
                meta.type_name
            },
            subtype: if meta.subtype.is_empty() {
                "AMR".to_string()
            } else {
                meta.subtype
            },
            reportable: meta.reportable,
            db_version: version.clone(),
            seq: record.seq,
        });
    }

    anyhow::ensure!(!references.is_empty(), "no AMR CDS references loaded");
    Ok(references)
}

fn fallback_element_symbol(
    meta_element_symbol: &str,
    gene_symbol: &str,
    allele_symbol: &str,
) -> String {
    if !meta_element_symbol.is_empty() {
        return meta_element_symbol.to_string();
    }
    if !allele_symbol.is_empty() {
        return allele_symbol.to_string();
    }
    gene_symbol.to_string()
}

fn parse_cds_header(record: &FastaRecord) -> CdsHeader {
    let pipe: Vec<&str> = record.id.split('|').collect();
    let desc = record.description.as_str();
    let desc_pipe: Vec<&str> = desc.split('|').collect();
    let fields = if desc_pipe.len() > pipe.len() {
        desc_pipe
    } else {
        pipe
    };

    let product = fields
        .get(6)
        .copied()
        .unwrap_or("")
        .split_whitespace()
        .next()
        .unwrap_or("")
        .replace('_', " ");

    CdsHeader {
        protein_accession: fields.first().copied().unwrap_or(&record.id).to_string(),
        nucleotide_accession: fields.get(1).copied().unwrap_or("").to_string(),
        gene_symbol: fields.get(4).copied().unwrap_or(&record.id).to_string(),
        allele_symbol: fields.get(5).copied().unwrap_or("").to_string(),
        product,
    }
}

fn fallback_family(meta_family: &str, gene_symbol: &str, allele_symbol: &str) -> String {
    if !meta_family.is_empty() {
        return meta_family.to_string();
    }
    let symbol = if allele_symbol.is_empty() {
        gene_symbol
    } else {
        allele_symbol
    };
    symbol
        .split(['_', '-'])
        .next()
        .filter(|s| !s.is_empty())
        .unwrap_or(symbol)
        .to_string()
}

fn load_metadata(db_dir: &Path) -> anyhow::Result<HashMap<String, Metadata>> {
    let mut map = HashMap::new();
    load_catalog_metadata(db_dir, &mut map)?;
    load_hierarchy_metadata(db_dir, &mut map)?;
    load_fam_metadata(db_dir, &mut map)?;
    Ok(map)
}

fn load_catalog_metadata(db_dir: &Path, map: &mut HashMap<String, Metadata>) -> anyhow::Result<()> {
    let path = db_dir.join("ReferenceGeneCatalog.txt");
    if !path.exists() {
        return Ok(());
    }
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(());
    };
    let columns = columns(header);
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let meta = Metadata {
            element_symbol: first_non_empty([
                field(&fields, &columns, "allele"),
                field(&fields, &columns, "gene_family"),
                field(&fields, &columns, "hierarchy_node"),
            ]),
            class_name: field(&fields, &columns, "class").to_string(),
            subclass: field(&fields, &columns, "subclass").to_string(),
            hierarchy_node: field(&fields, &columns, "hierarchy_node").to_string(),
            scope: field(&fields, &columns, "scope").to_string(),
            type_name: field(&fields, &columns, "type").to_string(),
            subtype: field(&fields, &columns, "subtype").to_string(),
            reportable: 0,
        };
        insert_if_present(map, field(&fields, &columns, "allele"), &meta);
        insert_if_present(map, field(&fields, &columns, "gene_family"), &meta);
        insert_if_present(
            map,
            field(&fields, &columns, "refseq_protein_accession"),
            &meta,
        );
        insert_if_present(
            map,
            field(&fields, &columns, "genbank_protein_accession"),
            &meta,
        );
        insert_if_present(map, field(&fields, &columns, "hierarchy_node"), &meta);
    }
    Ok(())
}

fn load_hierarchy_metadata(
    db_dir: &Path,
    map: &mut HashMap<String, Metadata>,
) -> anyhow::Result<()> {
    let path = db_dir.join("ReferenceGeneHierarchy.txt");
    if !path.exists() {
        return Ok(());
    }
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(());
    };
    let columns = columns(header);
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let meta = Metadata {
            element_symbol: first_non_empty([
                field(&fields, &columns, "symbol"),
                field(&fields, &columns, "node_id"),
                field(&fields, &columns, "name"),
            ]),
            class_name: field(&fields, &columns, "class").to_string(),
            subclass: field(&fields, &columns, "subclass").to_string(),
            hierarchy_node: field(&fields, &columns, "node_id").to_string(),
            scope: field(&fields, &columns, "scope").to_string(),
            type_name: field(&fields, &columns, "type").to_string(),
            subtype: field(&fields, &columns, "subtype").to_string(),
            reportable: 0,
        };
        insert_if_present(map, field(&fields, &columns, "node_id"), &meta);
        insert_if_present(map, field(&fields, &columns, "symbol"), &meta);
        insert_if_present(map, field(&fields, &columns, "prot_acc"), &meta);
    }
    Ok(())
}

fn load_fam_metadata(db_dir: &Path, map: &mut HashMap<String, Metadata>) -> anyhow::Result<()> {
    let path = db_dir.join("fam.tsv");
    if !path.exists() {
        return Ok(());
    }
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(());
    };
    let columns = columns(header);
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let meta = Metadata {
            element_symbol: first_non_empty([
                field(&fields, &columns, "gene_symbol"),
                field(&fields, &columns, "node_id"),
                field(&fields, &columns, "family_name"),
            ]),
            class_name: field(&fields, &columns, "class").to_string(),
            subclass: field(&fields, &columns, "subclass").to_string(),
            hierarchy_node: field(&fields, &columns, "node_id").to_string(),
            scope: String::new(),
            type_name: field(&fields, &columns, "type").to_string(),
            subtype: field(&fields, &columns, "subtype").to_string(),
            reportable: field(&fields, &columns, "reportable").parse().unwrap_or(0),
        };
        insert_if_present(map, field(&fields, &columns, "node_id"), &meta);
        insert_if_present(map, field(&fields, &columns, "gene_symbol"), &meta);
        insert_if_present(map, field(&fields, &columns, "hmm_id"), &meta);
    }
    Ok(())
}

fn columns(header: &str) -> HashMap<String, usize> {
    header
        .trim_start_matches('#')
        .split('\t')
        .enumerate()
        .map(|(idx, col)| (col.to_string(), idx))
        .collect()
}

fn field<'a>(fields: &'a [&str], columns: &HashMap<String, usize>, name: &str) -> &'a str {
    columns
        .get(name)
        .and_then(|idx| fields.get(*idx))
        .copied()
        .unwrap_or("")
        .trim()
}

fn first_non_empty(values: [&str; 3]) -> String {
    values
        .into_iter()
        .find(|value| !value.is_empty())
        .unwrap_or("")
        .to_string()
}

fn insert_if_present(map: &mut HashMap<String, Metadata>, key: &str, meta: &Metadata) {
    if !key.is_empty() {
        map.entry(key.to_string()).or_insert_with(|| meta.clone());
    }
}
