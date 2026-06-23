use crate::fasta::{FastaRecord, read_fasta};
use crate::translate::{DEFAULT_BACTERIAL_TRANSLATION_TABLE, translate_cds};
use anyhow::{Context, anyhow, ensure};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReferenceType {
    Amr,
    Stress,
    Virulence,
}

impl ReferenceType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Amr => "AMR",
            Self::Stress => "STRESS",
            Self::Virulence => "VIRULENCE",
        }
    }
}

/// Struct to reconstruct/represent the hierarchy of genes/proteins that AMRFinderPlus uses,
/// see: https://github.com/ncbi/amr/wiki/AMRFinderPlus-database#referencegenehierarchytxt and
/// https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy/ .
#[derive(Debug, Clone)]
pub struct HierarchyNode {
    pub node_id: String,
    pub parent_node_id: String,
    pub symbol: String,
    pub class_name: String,
    pub subclass: String,
    pub scope: String,
    pub type_name: String,
    pub subtype: String,
    pub reportable: u8,
}

/// Struct for representing the entries from the reference AMRFinderPlus database. See
/// https://github.com/ncbi/amr/wiki/AMRFinderPlus-database#referencegenecatalogtxt and
/// https://www.ncbi.nlm.nih.gov/pathogens/isolates#/refgene/ .
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
    pub hierarchy_path: Vec<HierarchyNode>,
    pub db_version: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct NodeMetadata {
    node_id: String,
    parent_node_id: String,
    symbol: String,
    class_name: String,
    subclass: String,
    scope: String,
    type_name: String,
    subtype: String,
    reportable: Option<u8>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct CatalogEntry {
    hierarchy_node: String,
    gene_family: String,
    allele_symbol: String,
    product_name: String,
    scope: String,
    type_name: String,
    subtype: String,
    class_name: String,
    subclass: String,
    refseq_nucleotide_accessions: BTreeSet<String>,
    genbank_nucleotide_accessions: BTreeSet<String>,
    genbank_protein_accessions: BTreeSet<String>,
}

#[derive(Debug, Clone, Default)]
struct CatalogRow<'a> {
    refseq_protein_accession: &'a str,
    genbank_protein_accession: &'a str,
    refseq_nucleotide_accession: &'a str,
    genbank_nucleotide_accession: &'a str,
    hierarchy_node: &'a str,
    gene_family: &'a str,
    allele_symbol: &'a str,
    product_name: &'a str,
    scope: &'a str,
    type_name: &'a str,
    subtype: &'a str,
    class_name: &'a str,
    subclass: &'a str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FastaSource {
    Cds,
    Protein,
}

impl FastaSource {
    fn label(self) -> &'static str {
        match self {
            Self::Cds => "CDS FASTA",
            Self::Protein => "protein FASTA",
        }
    }

    fn min_fields(self) -> usize {
        match self {
            Self::Cds => 7,
            Self::Protein => 10,
        }
    }

    fn node_index(self) -> usize {
        match self {
            Self::Cds => 4,
            Self::Protein => 3,
        }
    }

    fn parent_node_index(self) -> usize {
        match self {
            Self::Cds => 5,
            Self::Protein => 4,
        }
    }

    fn mechanism(self, fields: &[&str]) -> String {
        match self {
            Self::Cds => String::new(),
            Self::Protein => fields[5].to_string(),
        }
    }

    fn nucleotide_accession(self, fields: &[&str]) -> String {
        match self {
            Self::Cds => fields[1].to_string(),
            Self::Protein => String::new(),
        }
    }

    fn fusion_part(self, fields: &[&str]) -> anyhow::Result<usize> {
        let idx = match self {
            Self::Cds => 2,
            Self::Protein => 1,
        };
        parse_header_usize(fields[idx], self.label(), "fusion part")
    }

    fn total_fusion_parts(self, fields: &[&str]) -> anyhow::Result<usize> {
        let idx = match self {
            Self::Cds => 3,
            Self::Protein => 2,
        };
        parse_header_usize(fields[idx], self.label(), "total fusion parts")
    }

    fn product_field(self, fields: &[&str]) -> String {
        match self {
            Self::Cds => fields[6],
            Self::Protein => fields[9],
        }
        .replace('_', " ")
    }
}

#[derive(Debug, Clone)]
struct ReferenceHeader {
    protein_accession: String,
    nucleotide_accession: String,
    fusion_part: usize,
    total_fusion_parts: usize,
    node_token: String,
    parent_node_token: String,
    mechanism: String,
    product: String,
}

// Main functions follows now

/// For loading the db info
pub fn load_amrfinder_references(
    db_dir: &Path,
    include_types: &[ReferenceType],
) -> anyhow::Result<Vec<AmrReference>> {
    load_amrfinder_references_from_fasta(
        db_dir,
        &db_dir.join("AMR_CDS.fa"),
        FastaSource::Cds,
        include_types,
    )
}

pub fn load_amrfinder_protein_references(
    db_dir: &Path,
    include_types: &[ReferenceType],
) -> anyhow::Result<Vec<AmrReference>> {
    let protein_path = db_dir.join("AMRProt.fa");
    if protein_path.exists() {
        return load_amrfinder_references_from_fasta(
            db_dir,
            &protein_path,
            FastaSource::Protein,
            include_types,
        );
    }

    let mut references = load_amrfinder_references(db_dir, include_types)?;
    for reference in &mut references {
        reference.seq = translate_cds(&reference.seq, DEFAULT_BACTERIAL_TRANSLATION_TABLE);
    }
    references.retain(|reference| !reference.seq.is_empty());
    anyhow::ensure!(
        !references.is_empty(),
        "no protein AMR references loaded from AMR_CDS.fa fallback"
    );
    Ok(references)
}

fn load_amrfinder_references_from_fasta(
    db_dir: &Path,
    fasta_path: &Path,
    source: FastaSource,
    include_types: &[ReferenceType],
) -> anyhow::Result<Vec<AmrReference>> {
    let version = fs::read_to_string(db_dir.join("version.txt"))
        .unwrap_or_else(|_| "unknown".to_string())
        .trim()
        .to_string();
    let hierarchy_nodes = load_hierarchy_nodes(db_dir)?;
    let node_metadata = load_node_metadata(db_dir)?;
    let catalog_entries = load_catalog_entries(db_dir)?;
    let protein_reportability = load_protein_reportability(db_dir)?;
    let records = read_fasta(fasta_path).with_context(|| {
        format!(
            "read AMRFinderPlus {} {}",
            source.label(),
            fasta_path.display()
        )
    })?;

    let mut references = Vec::new();
    for record in records {
        let header = parse_reference_header(&record, source)?;
        let Some(catalog) = catalog_entries.get(&header.protein_accession) else {
            if source == FastaSource::Protein && header.mechanism == "mutation" {
                continue;
            }
            return Err(anyhow!(
                "missing ReferenceGeneCatalog.txt exact-reference entry for protein accession {} in {}",
                header.protein_accession,
                fasta_path.display()
            ));
        };

        ensure!(
            !catalog.hierarchy_node.is_empty(),
            "empty hierarchy_node in ReferenceGeneCatalog.txt for protein accession {}",
            header.protein_accession
        );
        let effective_node = resolve_effective_node(&header, catalog, &node_metadata)?;

        let meta = node_metadata.get(&effective_node).ok_or_else(|| {
            anyhow!(
                "missing ReferenceGeneHierarchy.txt node '{}' for protein accession {}",
                effective_node,
                header.protein_accession
            )
        })?;

        if !include_types
            .iter()
            .any(|kind| kind.as_str() == meta.type_name)
        {
            continue;
        }

        validate_nucleotide_accession(&header, catalog, fasta_path)?;

        let reportable = protein_reportability
            .get(&header.protein_accession)
            .copied()
            .or(meta.reportable)
            .unwrap_or(0);
        let ref_hierarchy_path = hierarchy_path(&effective_node, &hierarchy_nodes, meta);
        let gene_symbol = if catalog.gene_family.is_empty() {
            meta.symbol.clone()
        } else {
            catalog.gene_family.clone()
        };
        let allele_symbol = catalog.allele_symbol.clone();
        references.push(AmrReference {
            protein_accession: header.protein_accession,
            nucleotide_accession: header.nucleotide_accession,
            element_symbol: fallback_element_symbol(
                &allele_symbol,
                &meta.symbol,
                &gene_symbol,
                &effective_node,
            ),
            gene_symbol: gene_symbol.clone(),
            allele_symbol: allele_symbol.clone(),
            product: if catalog.product_name.is_empty() {
                header.product
            } else {
                catalog.product_name.clone()
            },
            family: fallback_family(&gene_symbol, &meta.symbol, &effective_node),
            class_name: meta.class_name.clone(),
            subclass: meta.subclass.clone(),
            hierarchy_node: effective_node,
            scope: meta.scope.clone(),
            type_name: meta.type_name.clone(),
            subtype: meta.subtype.clone(),
            reportable,
            hierarchy_path: ref_hierarchy_path,
            db_version: version.clone(),
            seq: record.seq,
        });
    }

    anyhow::ensure!(
        !references.is_empty(),
        "no selected AMRFinderPlus references loaded from {}",
        fasta_path.display()
    );
    Ok(references)
}


// Helper/other functions that serve to parse essentially, solve some issues that can happen, etc.
fn fallback_element_symbol(
    allele_symbol: &str,
    hierarchy_symbol: &str,
    gene_symbol: &str,
    hierarchy_node: &str,
) -> String {
    if !allele_symbol.is_empty() {
        return allele_symbol.to_string();
    }
    if !hierarchy_symbol.is_empty() {
        return hierarchy_symbol.to_string();
    }
    if !gene_symbol.is_empty() {
        return gene_symbol.to_string();
    }
    hierarchy_node.to_string()
}


fn fallback_family(gene_symbol: &str, hierarchy_symbol: &str, hierarchy_node: &str) -> String {
    if !gene_symbol.is_empty() {
        return gene_symbol.to_string();
    }
    if !hierarchy_symbol.is_empty() {
        return hierarchy_symbol.to_string();
    }
    hierarchy_node.to_string()
}


fn parse_reference_header(
    record: &FastaRecord,
    source: FastaSource,
) -> anyhow::Result<ReferenceHeader> {
    let fields: Vec<&str> = record.id.split('|').collect();
    ensure!(
        fields.len() >= source.min_fields(),
        "invalid {} header for {}: expected at least {} pipe-delimited fields, got {}",
        source.label(),
        record.id,
        source.min_fields(),
        fields.len()
    );

    Ok(ReferenceHeader {
        protein_accession: fields[0].to_string(),
        nucleotide_accession: source.nucleotide_accession(&fields),
        fusion_part: source.fusion_part(&fields)?,
        total_fusion_parts: source.total_fusion_parts(&fields)?,
        node_token: fields[source.node_index()].to_string(),
        parent_node_token: fields[source.parent_node_index()].to_string(),
        mechanism: source.mechanism(&fields),
        product: source.product_field(&fields),
    })
}

fn parse_header_usize(value: &str, source_label: &str, field_name: &str) -> anyhow::Result<usize> {
    value
        .parse::<usize>()
        .with_context(|| format!("parse {} {} value '{}'", source_label, field_name, value))
}


fn resolve_effective_node(
    header: &ReferenceHeader,
    catalog: &CatalogEntry,
    node_metadata: &HashMap<String, NodeMetadata>,
) -> anyhow::Result<String> {
    if header.node_token == catalog.hierarchy_node
        || header.parent_node_token == catalog.hierarchy_node
    {
        return Ok(catalog.hierarchy_node.clone());
    }

    let fused_parts: Vec<&str> = catalog.hierarchy_node.split(',').map(str::trim).collect();
    let header_matches_fused = fused_parts.iter().any(|part| {
        !part.is_empty() && (*part == header.node_token || *part == header.parent_node_token)
    });

    if header.total_fusion_parts > 1 && header_matches_fused {
        let candidate = if node_metadata.contains_key(&header.parent_node_token) {
            header.parent_node_token.clone()
        } else {
            header.node_token.clone()
        };
        if node_metadata.contains_key(&candidate) {
            return Ok(candidate);
        }
    }

    Err(anyhow!(
        "header/catalog node mismatch for protein accession {}: header node='{}', header parent='{}', catalog hierarchy_node='{}', fusion_part={}/{}",
        header.protein_accession,
        header.node_token,
        header.parent_node_token,
        catalog.hierarchy_node,
        header.fusion_part,
        header.total_fusion_parts
    ))
}


/// This just checks that there is no issue with a particular nt accession, by reviewing that there is an entry in the catalog
fn validate_nucleotide_accession(
    header: &ReferenceHeader,
    catalog: &CatalogEntry,
    fasta_path: &Path,
) -> anyhow::Result<()> {
    if header.nucleotide_accession.is_empty() {
        return Ok(());
    }
    if catalog
        .refseq_nucleotide_accessions
        .contains(&header.nucleotide_accession)
        || catalog
            .genbank_nucleotide_accessions
            .contains(&header.nucleotide_accession)
    {
        return Ok(());
    }
    Err(anyhow!(
        "header/catalog nucleotide accession mismatch for protein accession {} in {}: header='{}', refseq={:?}, genbank={:?}",
        header.protein_accession,
        fasta_path.display(),
        header.nucleotide_accession,
        catalog.refseq_nucleotide_accessions,
        catalog.genbank_nucleotide_accessions
    ))
}


fn load_node_metadata(db_dir: &Path) -> anyhow::Result<HashMap<String, NodeMetadata>> {
    let path = db_dir.join("ReferenceGeneHierarchy.txt");
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(HashMap::new());
    };
    let columns = columns(header);
    let reportability = load_fam_reportability(db_dir)?;
    let mut map = HashMap::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let node_id = field(&fields, &columns, "node_id");
        if node_id.is_empty() {
            continue;
        }
        map.insert(
            node_id.to_string(),
            NodeMetadata {
                node_id: node_id.to_string(),
                parent_node_id: field(&fields, &columns, "parent_node_id").to_string(),
                symbol: field(&fields, &columns, "symbol").to_string(),
                class_name: field(&fields, &columns, "class").to_string(),
                subclass: field(&fields, &columns, "subclass").to_string(),
                scope: field(&fields, &columns, "scope").to_string(),
                type_name: field(&fields, &columns, "type").to_string(),
                subtype: field(&fields, &columns, "subtype").to_string(),
                reportable: reportability.get(node_id).copied(),
            },
        );
    }
    Ok(map)
}


fn load_fam_reportability(db_dir: &Path) -> anyhow::Result<HashMap<String, u8>> {
    let path = db_dir.join("fam.tsv");
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(HashMap::new());
    };
    let columns = columns(header);
    let mut reportability = HashMap::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let node_id = field(&fields, &columns, "node_id");
        if node_id.is_empty() {
            continue;
        }
        let value = field(&fields, &columns, "reportable");
        if value.is_empty() {
            continue;
        }
        let parsed = value.parse::<u8>().with_context(|| {
            format!(
                "parse fam.tsv reportable value '{}' for node {}",
                value, node_id
            )
        })?;
        reportability.insert(node_id.to_string(), parsed);
    }
    Ok(reportability)
}


fn load_protein_reportability(db_dir: &Path) -> anyhow::Result<HashMap<String, u8>> {
    let path = db_dir.join("AMRProt.fa");
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let records = read_fasta(&path).with_context(|| format!("read FASTA {}", path.display()))?;
    let mut reportability = HashMap::new();
    for record in records {
        let fields: Vec<&str> = record.id.split('|').collect();
        if fields.len() < FastaSource::Protein.min_fields() {
            continue;
        }
        let protein_accession = fields[0].trim();
        let value = fields[6].trim();
        if protein_accession.is_empty() || value.is_empty() {
            continue;
        }
        let parsed = value.parse::<u8>().with_context(|| {
            format!(
                "parse AMRProt.fa reportable value '{}' for protein {}",
                value, protein_accession
            )
        })?;
        if let Some(existing) = reportability.insert(protein_accession.to_string(), parsed)
            && existing != parsed
        {
            return Err(anyhow!(
                "conflicting AMRProt.fa reportable values for protein accession {}: {} vs {}",
                protein_accession,
                existing,
                parsed
            ));
        }
    }
    Ok(reportability)
}


fn load_catalog_entries(db_dir: &Path) -> anyhow::Result<HashMap<String, CatalogEntry>> {
    let path = db_dir.join("ReferenceGeneCatalog.txt");
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(HashMap::new());
    };
    let columns = columns(header);
    let mut entries_by_refseq: HashMap<String, CatalogEntry> = HashMap::new();
    let mut entries_by_genbank: HashMap<String, CatalogEntry> = HashMap::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let row = CatalogRow {
            refseq_protein_accession: field(&fields, &columns, "refseq_protein_accession"),
            genbank_protein_accession: field(&fields, &columns, "genbank_protein_accession"),
            refseq_nucleotide_accession: field(&fields, &columns, "refseq_nucleotide_accession"),
            genbank_nucleotide_accession: field(&fields, &columns, "genbank_nucleotide_accession"),
            hierarchy_node: field(&fields, &columns, "hierarchy_node"),
            gene_family: field(&fields, &columns, "gene_family"),
            allele_symbol: field(&fields, &columns, "allele"),
            product_name: field(&fields, &columns, "product_name"),
            scope: field(&fields, &columns, "scope"),
            type_name: field(&fields, &columns, "type"),
            subtype: field(&fields, &columns, "subtype"),
            class_name: field(&fields, &columns, "class"),
            subclass: field(&fields, &columns, "subclass"),
        };
        if !is_exact_reference_catalog_row(&row) {
            continue;
        }
        let entry = CatalogEntry::from_row(&row);
        insert_catalog_entry(&mut entries_by_refseq, row.refseq_protein_accession, &entry)?;
        insert_catalog_entry(
            &mut entries_by_genbank,
            row.genbank_protein_accession,
            &entry,
        )?;
    }

    let mut map = HashMap::new();
    for (key, entry) in entries_by_refseq.into_iter().chain(entries_by_genbank) {
        insert_catalog_alias(&mut map, &key, &entry)?;
    }
    Ok(map)
}

fn load_hierarchy_nodes(db_dir: &Path) -> anyhow::Result<HashMap<String, HierarchyNode>> {
    let path = db_dir.join("ReferenceGeneHierarchy.txt");
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let text = fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(HashMap::new());
    };
    let columns = columns(header);
    let reportability = load_fam_reportability(db_dir)?;
    let mut nodes = HashMap::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let node_id = field(&fields, &columns, "node_id");
        if node_id.is_empty() {
            continue;
        }

        nodes.insert(
            node_id.to_string(),
            HierarchyNode {
                node_id: node_id.to_string(),
                parent_node_id: field(&fields, &columns, "parent_node_id").to_string(),
                symbol: first_non_empty([
                    field(&fields, &columns, "symbol"),
                    node_id,
                    field(&fields, &columns, "name"),
                ]),
                class_name: field(&fields, &columns, "class").to_string(),
                subclass: field(&fields, &columns, "subclass").to_string(),
                scope: field(&fields, &columns, "scope").to_string(),
                type_name: field(&fields, &columns, "type").to_string(),
                subtype: field(&fields, &columns, "subtype").to_string(),
                reportable: reportability.get(node_id).copied().unwrap_or(0),
            },
        );
    }
    Ok(nodes)
}


/// This gives you the branch as a vector of nodes up to the top
fn hierarchy_path(
    leaf_node_id: &str,
    nodes: &HashMap<String, HierarchyNode>,
    meta: &NodeMetadata,
) -> Vec<HierarchyNode> {
    let mut path = Vec::new();
    let mut seen = HashSet::new();
    let mut current = leaf_node_id;
    while !current.is_empty() && seen.insert(current.to_string()) {
        let Some(node) = nodes.get(current) else {
            break;
        };
        path.push(node.clone());
        current = &node.parent_node_id;
    }
    if path.is_empty() && !leaf_node_id.is_empty() {
        path.push(HierarchyNode {
            node_id: leaf_node_id.to_string(),
            parent_node_id: String::new(),
            symbol: first_non_empty([&meta.symbol, leaf_node_id, ""]),
            class_name: meta.class_name.clone(),
            subclass: meta.subclass.clone(),
            scope: meta.scope.clone(),
            type_name: meta.type_name.clone(),
            subtype: meta.subtype.clone(),
            reportable: meta.reportable.unwrap_or(0),
        });
    }
    path
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

impl CatalogEntry {
    fn from_row(row: &CatalogRow<'_>) -> Self {
        let mut entry = Self {
            hierarchy_node: row.hierarchy_node.to_string(),
            gene_family: row.gene_family.to_string(),
            allele_symbol: row.allele_symbol.to_string(),
            product_name: row.product_name.to_string(),
            scope: row.scope.to_string(),
            type_name: row.type_name.to_string(),
            subtype: row.subtype.to_string(),
            class_name: row.class_name.to_string(),
            subclass: row.subclass.to_string(),
            refseq_nucleotide_accessions: BTreeSet::new(),
            genbank_nucleotide_accessions: BTreeSet::new(),
            genbank_protein_accessions: BTreeSet::new(),
        };
        if !row.refseq_nucleotide_accession.is_empty() {
            entry
                .refseq_nucleotide_accessions
                .insert(row.refseq_nucleotide_accession.to_string());
        }
        if !row.genbank_nucleotide_accession.is_empty() {
            entry
                .genbank_nucleotide_accessions
                .insert(row.genbank_nucleotide_accession.to_string());
        }
        if !row.genbank_protein_accession.is_empty() {
            entry
                .genbank_protein_accessions
                .insert(row.genbank_protein_accession.to_string());
        }
        entry
    }
}

fn is_exact_reference_catalog_row(row: &CatalogRow<'_>) -> bool {
    !row.hierarchy_node.is_empty() && !row.subtype.starts_with("POINT")
}

fn insert_catalog_entry(
    map: &mut HashMap<String, CatalogEntry>,
    key: &str,
    entry: &CatalogEntry,
) -> anyhow::Result<()> {
    if key.is_empty() {
        return Ok(());
    }
    if let Some(existing) = map.get_mut(key) {
        merge_catalog_entry(existing, entry).with_context(|| {
            format!(
                "conflicting ReferenceGeneCatalog.txt entries for protein accession {}",
                key
            )
        })?;
        return Ok(());
    }
    map.insert(key.to_string(), entry.clone());
    Ok(())
}

fn merge_catalog_entry(existing: &mut CatalogEntry, incoming: &CatalogEntry) -> anyhow::Result<()> {
    if existing.hierarchy_node != incoming.hierarchy_node
        || existing.gene_family != incoming.gene_family
        || existing.allele_symbol != incoming.allele_symbol
        || existing.product_name != incoming.product_name
        || existing.scope != incoming.scope
        || existing.type_name != incoming.type_name
        || existing.subtype != incoming.subtype
        || existing.class_name != incoming.class_name
        || existing.subclass != incoming.subclass
    {
        return Err(anyhow!("existing={:?}, incoming={:?}", existing, incoming));
    }
    existing
        .refseq_nucleotide_accessions
        .extend(incoming.refseq_nucleotide_accessions.iter().cloned());
    existing
        .genbank_nucleotide_accessions
        .extend(incoming.genbank_nucleotide_accessions.iter().cloned());
    existing
        .genbank_protein_accessions
        .extend(incoming.genbank_protein_accessions.iter().cloned());
    Ok(())
}

fn insert_catalog_alias(
    map: &mut HashMap<String, CatalogEntry>,
    key: &str,
    entry: &CatalogEntry,
) -> anyhow::Result<()> {
    if key.is_empty() {
        return Ok(());
    }
    if let Some(existing) = map.get(key)
        && existing != entry
    {
        return Err(anyhow!(
            "conflicting merged catalog aliases for protein accession {}: existing={:?}, incoming={:?}",
            key,
            existing,
            entry
        ));
    }
    map.insert(key.to_string(), entry.clone());
    Ok(())
}



// ======================================================== TESTS
#[cfg(test)]
mod tests {
    use super::*;

    fn sample_catalog_entry() -> CatalogEntry {
        let mut genbank_nucleotide_accessions = BTreeSet::new();
        genbank_nucleotide_accessions.insert("L11078.1".to_string());
        CatalogEntry {
            hierarchy_node: "stxA2b".to_string(),
            gene_family: "stxA2".to_string(),
            allele_symbol: String::new(),
            product_name: "Shiga toxin Stx2b subunit A".to_string(),
            scope: "plus".to_string(),
            type_name: "VIRULENCE".to_string(),
            subtype: "VIRULENCE".to_string(),
            class_name: "STX2".to_string(),
            subclass: "stxA2".to_string(),
            refseq_nucleotide_accessions: BTreeSet::new(),
            genbank_nucleotide_accessions,
            genbank_protein_accessions: BTreeSet::new(),
        }
    }

    fn sample_node_metadata() -> NodeMetadata {
        NodeMetadata {
            node_id: "stxA2b".to_string(),
            parent_node_id: "stxA2".to_string(),
            symbol: "stxA2b".to_string(),
            class_name: "STX2".to_string(),
            subclass: "stxA2b".to_string(),
            scope: "plus".to_string(),
            type_name: "VIRULENCE".to_string(),
            subtype: "VIRULENCE".to_string(),
            reportable: Some(1),
        }
    }

    fn sample_catalog_row<'a>(refseq_protein_accession: &'a str) -> CatalogRow<'a> {
        CatalogRow {
            refseq_protein_accession,
            genbank_protein_accession: "AAS07596.1",
            refseq_nucleotide_accession: "",
            genbank_nucleotide_accession: "L11078.1",
            hierarchy_node: "stxA2b",
            gene_family: "stxA2",
            allele_symbol: "",
            product_name: "Shiga toxin Stx2b subunit A",
            scope: "plus",
            type_name: "VIRULENCE",
            subtype: "VIRULENCE",
            class_name: "STX2",
            subclass: "stxA2",
        }
    }

    #[test]
    fn parses_cds_header_layout() {
        let record = FastaRecord {
            id: "AAA16360.1|L11078.1|1|1|stxA2b|stxA2b|Shiga_toxin_Stx2b_subunit_A".to_string(),
            description: "AAA16360.1|L11078.1|1|1|stxA2b|stxA2b|Shiga_toxin_Stx2b_subunit_A L11078.1:177-1136"
                .to_string(),
            seq: Vec::new(),
        };
        let header = parse_reference_header(&record, FastaSource::Cds).unwrap();
        assert_eq!(header.protein_accession, "AAA16360.1");
        assert_eq!(header.nucleotide_accession, "L11078.1");
        assert_eq!(header.fusion_part, 1);
        assert_eq!(header.total_fusion_parts, 1);
        assert_eq!(header.node_token, "stxA2b");
        assert_eq!(header.parent_node_token, "stxA2b");
        assert_eq!(header.mechanism, "");
        assert_eq!(header.product, "Shiga toxin Stx2b subunit A");
    }

    #[test]
    fn parses_protein_header_layout() {
        let record = FastaRecord {
            id: "AAA16360.1|1|1|stxA2b|stxA2b||1|stxA2b|STX2|Shiga_toxin_Stx2b_subunit_A"
                .to_string(),
            description: "AAA16360.1|1|1|stxA2b|stxA2b||1|stxA2b|STX2|Shiga_toxin_Stx2b_subunit_A"
                .to_string(),
            seq: Vec::new(),
        };
        let header = parse_reference_header(&record, FastaSource::Protein).unwrap();
        assert_eq!(header.protein_accession, "AAA16360.1");
        assert_eq!(header.nucleotide_accession, "");
        assert_eq!(header.fusion_part, 1);
        assert_eq!(header.total_fusion_parts, 1);
        assert_eq!(header.node_token, "stxA2b");
        assert_eq!(header.parent_node_token, "stxA2b");
        assert_eq!(header.mechanism, "");
        assert_eq!(header.product, "Shiga toxin Stx2b subunit A");
    }

    #[test]
    fn parses_protein_fusion_part_layout() {
        let record = FastaRecord {
            id: "WP_001028144.1|1|2|aac(6')-Ie|aac(6')-Ie|acetyltransferase|2|AMIKACIN/KANAMYCIN/TOBRAMYCIN|AMINOGLYCOSIDE|bifunctional_aminoglycoside_N-acetyltransferase_AAC(6')-Ie/aminoglycoside_O-phosphotransferase_APH(2'')-Ia"
                .to_string(),
            description: String::new(),
            seq: Vec::new(),
        };
        let header = parse_reference_header(&record, FastaSource::Protein).unwrap();
        assert_eq!(header.protein_accession, "WP_001028144.1");
        assert_eq!(header.fusion_part, 1);
        assert_eq!(header.total_fusion_parts, 2);
        assert_eq!(header.node_token, "aac(6')-Ie");
        assert_eq!(header.parent_node_token, "aac(6')-Ie");
    }

    #[test]
    fn parses_protein_mutation_mechanism() {
        let record = FastaRecord {
            id: "WP_000019358.1|1|1|soxS|soxS|mutation|2|||regulatory_protein_SoxS".to_string(),
            description: String::new(),
            seq: Vec::new(),
        };
        let header = parse_reference_header(&record, FastaSource::Protein).unwrap();
        assert_eq!(header.protein_accession, "WP_000019358.1");
        assert_eq!(header.mechanism, "mutation");
        assert_eq!(header.node_token, "soxS");
    }

    #[test]
    fn fallback_element_symbol_prefers_allele_then_hierarchy_symbol() {
        assert_eq!(
            fallback_element_symbol("blaTEM-1", "blaTEM", "blaTEM", "blaTEM-1"),
            "blaTEM-1"
        );
        assert_eq!(
            fallback_element_symbol("", "blaTEM", "", "blaTEM-1"),
            "blaTEM"
        );
    }

    #[test]
    fn catalog_alias_insert_accepts_identical_values() {
        let mut map = HashMap::new();
        let entry = sample_catalog_entry();
        insert_catalog_alias(&mut map, "AAA16360.1", &entry).unwrap();
        insert_catalog_alias(&mut map, "AAA16360.1", &entry).unwrap();
        assert_eq!(map.get("AAA16360.1"), Some(&entry));
    }

    #[test]
    fn catalog_alias_insert_rejects_conflicting_values() {
        let mut map = HashMap::new();
        let mut entry = sample_catalog_entry();
        insert_catalog_alias(&mut map, "AAA16360.1", &entry).unwrap();
        entry.hierarchy_node = "stxA2".to_string();
        let err = insert_catalog_alias(&mut map, "AAA16360.1", &entry).unwrap_err();
        assert!(
            err.to_string()
                .contains("conflicting merged catalog aliases")
        );
    }

    #[test]
    fn validates_matching_nucleotide_accession() {
        let header = ReferenceHeader {
            protein_accession: "AAA16360.1".to_string(),
            nucleotide_accession: "L11078.1".to_string(),
            fusion_part: 1,
            total_fusion_parts: 1,
            node_token: "stxA2b".to_string(),
            parent_node_token: "stxA2b".to_string(),
            mechanism: String::new(),
            product: String::new(),
        };
        validate_nucleotide_accession(&header, &sample_catalog_entry(), Path::new("db")).unwrap();
    }

    #[test]
    fn reports_mismatching_nucleotide_accession() {
        let header = ReferenceHeader {
            protein_accession: "AAA16360.1".to_string(),
            nucleotide_accession: "WRONG".to_string(),
            fusion_part: 1,
            total_fusion_parts: 1,
            node_token: "stxA2b".to_string(),
            parent_node_token: "stxA2b".to_string(),
            mechanism: String::new(),
            product: String::new(),
        };
        let err = validate_nucleotide_accession(&header, &sample_catalog_entry(), Path::new("db"))
            .unwrap_err();
        assert!(
            err.to_string()
                .contains("header/catalog nucleotide accession mismatch")
        );
    }

    #[test]
    fn hierarchy_path_uses_node_metadata_fallback() {
        let path = hierarchy_path("stxA2b", &HashMap::new(), &sample_node_metadata());
        assert_eq!(path.len(), 1);
        assert_eq!(path[0].node_id, "stxA2b");
        assert_eq!(path[0].symbol, "stxA2b");
    }

    #[test]
    fn exact_catalog_rows_merge_when_only_accessions_differ() {
        let mut existing = CatalogEntry::from_row(&CatalogRow {
            genbank_nucleotide_accession: "AY443052.1",
            ..sample_catalog_row("WP_000649751.1")
        });
        let incoming = CatalogEntry::from_row(&CatalogRow {
            genbank_nucleotide_accession: "AB015057.1",
            genbank_protein_accession: "BAA34372.1",
            ..sample_catalog_row("WP_000649751.1")
        });
        merge_catalog_entry(&mut existing, &incoming).unwrap();
        assert!(
            existing
                .genbank_nucleotide_accessions
                .contains("AY443052.1")
        );
        assert!(
            existing
                .genbank_nucleotide_accessions
                .contains("AB015057.1")
        );
        assert!(existing.genbank_protein_accessions.contains("BAA34372.1"));
    }

    #[test]
    fn point_mutation_catalog_rows_are_excluded() {
        let exact = sample_catalog_row("WP_015585966.1");
        let point = CatalogRow {
            hierarchy_node: "",
            allele_symbol: "fexA_G33A",
            subtype: "POINT",
            ..sample_catalog_row("WP_015585966.1")
        };
        let point_disrupt = CatalogRow {
            hierarchy_node: "acrR",
            allele_symbol: "",
            subtype: "POINT_DISRUPT",
            ..sample_catalog_row("WP_000101737.1")
        };
        assert!(is_exact_reference_catalog_row(&exact));
        assert!(!is_exact_reference_catalog_row(&point));
        assert!(!is_exact_reference_catalog_row(&point_disrupt));
    }

    #[test]
    fn resolves_fusion_part_to_component_node() {
        let mut node_metadata = HashMap::new();
        node_metadata.insert(
            "aac(6')-Ie".to_string(),
            NodeMetadata {
                node_id: "aac(6')-Ie".to_string(),
                parent_node_id: "aac(6')-Ie_fam".to_string(),
                symbol: "aac(6')-Ie".to_string(),
                class_name: "AMINOGLYCOSIDE".to_string(),
                subclass: "AMIKACIN/KANAMYCIN/TOBRAMYCIN".to_string(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: Some(2),
            },
        );
        let header = ReferenceHeader {
            protein_accession: "WP_001028144.1".to_string(),
            nucleotide_accession: "NG_047055.1".to_string(),
            fusion_part: 1,
            total_fusion_parts: 2,
            node_token: "aac(6')-Ie".to_string(),
            parent_node_token: "aac(6')-Ie".to_string(),
            mechanism: String::new(),
            product: String::new(),
        };
        let mut catalog = sample_catalog_entry();
        catalog.hierarchy_node = "aac(6')-Ie,aph(2'')-Ia".to_string();
        let resolved = resolve_effective_node(&header, &catalog, &node_metadata).unwrap();
        assert_eq!(resolved, "aac(6')-Ie");
    }

    #[test]
    fn merge_catalog_entry_rejects_biological_conflicts() {
        let mut existing = CatalogEntry::from_row(&sample_catalog_row("WP_000649751.1"));
        let incoming = CatalogEntry::from_row(&CatalogRow {
            hierarchy_node: "different_node",
            ..sample_catalog_row("WP_000649751.1")
        });
        let err = merge_catalog_entry(&mut existing, &incoming).unwrap_err();
        assert!(err.to_string().contains("existing="));
    }
}
