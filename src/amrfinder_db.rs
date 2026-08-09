use crate::fasta::{FastaRecord, read_fasta};
use crate::translate::{DEFAULT_BACTERIAL_TRANSLATION_TABLE, translate_cds};
use anyhow::{Context, anyhow, ensure};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
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
    let mut fusion_groups: BTreeMap<String, Vec<(usize, AmrReference)>> = BTreeMap::new();
    for record in records {
        let header = parse_reference_header(&record, source)?;
        let fusion_part = header.fusion_part;
        let total_fusion_parts = header.total_fusion_parts;
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
        let ref_hierarchy_path = hierarchy_path(&effective_node, &hierarchy_nodes, meta)?;
        let gene_symbol = if catalog.gene_family.is_empty() {
            meta.symbol.clone()
        } else {
            catalog.gene_family.clone()
        };
        let allele_symbol = catalog.allele_symbol.clone();
        let reference = AmrReference {
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
        };
        if total_fusion_parts > 1 {
            fusion_groups
                .entry(reference.protein_accession.clone())
                .or_default()
                .push((fusion_part, reference));
        } else {
            references.push(reference);
        }
    }

    for (accession, mut parts) in fusion_groups {
        parts.sort_by_key(|(part, _)| *part);
        let identical_seqs = parts.windows(2).all(|pair| pair[0].1.seq == pair[1].1.seq);
        if !identical_seqs {
            // A future DB could ship genuinely part-sliced sequences; those stay separate genes.
            references.extend(parts.into_iter().map(|(_, reference)| reference));
            continue;
        }
        let catalog = catalog_entries.get(&accession);
        references.push(merge_fusion_parts(parts, catalog)?);
    }

    anyhow::ensure!(
        !references.is_empty(),
        "no selected AMRFinderPlus references loaded from {}",
        fasta_path.display()
    );
    Ok(references)
}

/// AMRFinderPlus ships a fusion gene as one FASTA record per fused part, each tagged
/// `fusion_part|total_fusion_parts` in the header but carrying the full-length sequence.
/// Collapse identical-sequence parts into a single reference: kept separate, the twins make
/// every k-mer ambiguous and the fusion can never be reported.
fn merge_fusion_parts(
    parts: Vec<(usize, AmrReference)>,
    catalog: Option<&CatalogEntry>,
) -> anyhow::Result<AmrReference> {
    ensure!(!parts.is_empty(), "empty fusion part group");
    let parts: Vec<AmrReference> = parts.into_iter().map(|(_, reference)| reference).collect();
    if parts.len() == 1 {
        return Ok(parts.into_iter().next().expect("single fusion part"));
    }

    let element_symbol = join_distinct(parts.iter().map(|part| part.element_symbol.as_str()));
    let joined_gene_symbol = join_distinct(parts.iter().map(|part| part.gene_symbol.as_str()));
    let hierarchy_path = merge_fusion_paths(&parts);
    let reportable = parts.iter().map(|part| part.reportable).max().unwrap_or(0);
    let first = &parts[0];

    let catalog_or = |value: Option<&str>, fallback: &str| -> String {
        match value {
            Some(text) if !text.is_empty() => text.to_string(),
            _ => fallback.to_string(),
        }
    };

    Ok(AmrReference {
        protein_accession: first.protein_accession.clone(),
        nucleotide_accession: first.nucleotide_accession.clone(),
        element_symbol,
        gene_symbol: catalog_or(
            catalog.map(|entry| entry.gene_family.as_str()),
            &joined_gene_symbol,
        ),
        allele_symbol: catalog
            .map(|entry| entry.allele_symbol.clone())
            .unwrap_or_default(),
        product: first.product.clone(),
        family: catalog_or(
            catalog.map(|entry| entry.gene_family.as_str()),
            &first.family,
        ),
        class_name: catalog_or(
            catalog.map(|entry| entry.class_name.as_str()),
            &first.class_name,
        ),
        subclass: catalog_or(
            catalog.map(|entry| entry.subclass.as_str()),
            &first.subclass,
        ),
        hierarchy_node: catalog_or(
            catalog.map(|entry| entry.hierarchy_node.as_str()),
            &first.hierarchy_node,
        ),
        scope: first.scope.clone(),
        type_name: catalog_or(
            catalog.map(|entry| entry.type_name.as_str()),
            &first.type_name,
        ),
        subtype: catalog_or(catalog.map(|entry| entry.subtype.as_str()), &first.subtype),
        reportable,
        hierarchy_path,
        db_version: first.db_version.clone(),
        seq: first.seq.clone(),
    })
}

fn join_distinct<'a>(values: impl Iterator<Item = &'a str>) -> String {
    let mut joined: Vec<&str> = Vec::new();
    for value in values {
        if !value.is_empty() && !joined.contains(&value) {
            joined.push(value);
        }
    }
    joined.join("/")
}

/// Merged fusion path: each part's own nodes (leaf-first, per part order), then the ancestors
/// shared by every part, once, in the first part's order. Keeping leaf-first ordering lets
/// k-mers shared with either family's relatives resolve to that family's hierarchy units.
fn merge_fusion_paths(parts: &[AmrReference]) -> Vec<HierarchyNode> {
    let shared: Vec<HierarchyNode> = parts[0]
        .hierarchy_path
        .iter()
        .filter(|node| {
            parts[1..].iter().all(|part| {
                part.hierarchy_path
                    .iter()
                    .any(|other| other.node_id == node.node_id)
            })
        })
        .cloned()
        .collect();
    let shared_ids: HashSet<&str> = shared.iter().map(|node| node.node_id.as_str()).collect();

    let mut merged = Vec::new();
    let mut seen = HashSet::new();
    for part in parts {
        for node in &part.hierarchy_path {
            if shared_ids.contains(node.node_id.as_str()) {
                continue;
            }
            if seen.insert(node.node_id.clone()) {
                merged.push(node.clone());
            }
        }
    }
    merged.extend(shared);
    merged
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
        let type_name = field(&fields, &columns, "type");
        let subtype = field(&fields, &columns, "subtype");
        map.insert(
            node_id.to_string(),
            NodeMetadata {
                node_id: node_id.to_string(),
                parent_node_id: field(&fields, &columns, "parent_node_id").to_string(),
                symbol: field(&fields, &columns, "symbol").to_string(),
                class_name: plus_label_fallback(
                    type_name,
                    field(&fields, &columns, "class"),
                    subtype,
                ),
                subclass: plus_label_fallback(
                    type_name,
                    field(&fields, &columns, "subclass"),
                    subtype,
                ),
                scope: field(&fields, &columns, "scope").to_string(),
                type_name: type_name.to_string(),
                subtype: subtype.to_string(),
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
        let entry = CatalogEntry::from_row(&row)?;
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

        let type_name = field(&fields, &columns, "type");
        let subtype = field(&fields, &columns, "subtype");
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
                class_name: plus_label_fallback(
                    type_name,
                    field(&fields, &columns, "class"),
                    subtype,
                ),
                subclass: plus_label_fallback(
                    type_name,
                    field(&fields, &columns, "subclass"),
                    subtype,
                ),
                scope: field(&fields, &columns, "scope").to_string(),
                type_name: type_name.to_string(),
                subtype: subtype.to_string(),
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
) -> anyhow::Result<Vec<HierarchyNode>> {
    ensure!(
        !leaf_node_id.trim().is_empty(),
        "empty hierarchy leaf node for metadata {:?}",
        meta
    );

    let mut path = Vec::new();
    let mut seen = HashSet::new();
    let mut current = leaf_node_id;
    while !current.is_empty() {
        ensure!(
            seen.insert(current.to_string()),
            "cycle in ReferenceGeneHierarchy.txt while walking '{}': {}",
            leaf_node_id,
            hierarchy_diagnostic(current, nodes)
        );
        let Some(node) = nodes.get(current) else {
            return Err(anyhow!(
                "missing hierarchy node '{}' while walking '{}': {}",
                current,
                leaf_node_id,
                hierarchy_diagnostic(current, nodes)
            ));
        };
        path.push(node.clone());
        current = &node.parent_node_id;
    }

    ensure!(
        !path.is_empty(),
        "empty hierarchy path for '{}': {}",
        leaf_node_id,
        hierarchy_diagnostic(leaf_node_id, nodes)
    );
    Ok(path)
}

fn hierarchy_diagnostic(node_id: &str, nodes: &HashMap<String, HierarchyNode>) -> String {
    let mut children: HashMap<&str, Vec<&HierarchyNode>> = HashMap::new();
    for node in nodes.values() {
        children
            .entry(node.parent_node_id.as_str())
            .or_default()
            .push(node);
    }

    let node = nodes.get(node_id);
    let parent_id = node.map(|node| node.parent_node_id.as_str()).unwrap_or("");
    let parent = if parent_id.is_empty() {
        None
    } else {
        nodes.get(parent_id)
    };

    let mut siblings = if parent_id.is_empty() {
        Vec::new()
    } else {
        children.get(parent_id).cloned().unwrap_or_default()
    };
    siblings.sort_by(|a, b| a.node_id.cmp(&b.node_id));

    let mut descendants = Vec::new();
    let mut seen = HashSet::new();
    collect_descendants(node_id, &children, &mut descendants, &mut seen);
    descendants.sort_by(|a, b| a.node_id.cmp(&b.node_id));

    format!(
        "node={:?}; parent={:?}; siblings={:?}; descendants={:?}",
        node.map(format_node),
        parent.map(format_node),
        siblings.into_iter().map(format_node).collect::<Vec<_>>(),
        descendants.into_iter().map(format_node).collect::<Vec<_>>()
    )
}

fn collect_descendants<'a>(
    node_id: &str,
    children: &HashMap<&'a str, Vec<&'a HierarchyNode>>,
    out: &mut Vec<&'a HierarchyNode>,
    seen: &mut HashSet<String>,
) {
    if !seen.insert(node_id.to_string()) {
        return;
    }
    if let Some(kids) = children.get(node_id) {
        for child in kids {
            out.push(*child);
            collect_descendants(child.node_id.as_str(), children, out, seen);
        }
    }
}

fn format_node(node: &HierarchyNode) -> String {
    format!(
        "{{node_id='{}', parent='{}', symbol='{}', type='{}', subtype='{}', class='{}', subclass='{}', reportable={}}}",
        node.node_id,
        node.parent_node_id,
        node.symbol,
        node.type_name,
        node.subtype,
        node.class_name,
        node.subclass,
        node.reportable
    )
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

fn is_missing_amrfinder_label(value: &str) -> bool {
    let trimmed = value.trim();
    trimmed.is_empty()
        || trimmed.eq_ignore_ascii_case("NA")
        || trimmed.eq_ignore_ascii_case("Unclassified")
}

fn plus_label_fallback(type_name: &str, label: &str, subtype: &str) -> String {
    if is_missing_amrfinder_label(label)
        && !type_name.eq_ignore_ascii_case("AMR")
        && !is_missing_amrfinder_label(subtype)
    {
        subtype.trim().to_string()
    } else {
        label.trim().to_string()
    }
}

impl CatalogEntry {
    fn from_row(row: &CatalogRow<'_>) -> anyhow::Result<Self> {
        let class_name = plus_label_fallback(row.type_name, row.class_name, row.subtype);
        let subclass = plus_label_fallback(row.type_name, row.subclass, row.subtype);
        let mut entry = Self {
            hierarchy_node: row.hierarchy_node.to_string(),
            gene_family: row.gene_family.to_string(),
            allele_symbol: row.allele_symbol.to_string(),
            product_name: row.product_name.to_string(),
            scope: row.scope.to_string(),
            type_name: row.type_name.to_string(),
            subtype: row.subtype.to_string(),
            class_name,
            subclass,
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
        Ok(entry)
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
    fn hierarchy_path_rejects_missing_node_with_diagnostics() {
        let err = hierarchy_path("stxA2b", &HashMap::new(), &sample_node_metadata()).unwrap_err();
        let message = err.to_string();
        assert!(message.contains("missing hierarchy node 'stxA2b'"));
        assert!(message.contains("node="));
        assert!(message.contains("parent="));
        assert!(message.contains("siblings="));
        assert!(message.contains("descendants="));
    }

    #[test]
    fn exact_catalog_rows_merge_when_only_accessions_differ() {
        let mut existing = CatalogEntry::from_row(&CatalogRow {
            genbank_nucleotide_accession: "AY443052.1",
            ..sample_catalog_row("WP_000649751.1")
        })
        .unwrap();
        let incoming = CatalogEntry::from_row(&CatalogRow {
            genbank_nucleotide_accession: "AB015057.1",
            genbank_protein_accession: "BAA34372.1",
            ..sample_catalog_row("WP_000649751.1")
        })
        .unwrap();
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

    fn scratch_db_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "sparrowhawk_amr_dbtest_{}_{}",
            std::process::id(),
            tag
        ));
        if dir.exists() {
            fs::remove_dir_all(&dir).unwrap();
        }
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// Minimal on-disk DB with one two-part fusion gene (WP_F) and one ordinary
    /// relative (WP_R): partA -> famA -> ROOTX and partB -> famB -> ROOTX.
    fn write_fusion_test_db(dir: &Path, part1_seq: &str, part2_seq: &str) {
        fs::write(dir.join("version.txt"), "testdb\n").unwrap();
        fs::write(
            dir.join("ReferenceGeneHierarchy.txt"),
            "node_id\tparent_node_id\tsymbol\tclass\tsubclass\tscope\ttype\tsubtype\n\
             partA\tfamA\tpartA\tCLASSA\tSUBA\tcore\tAMR\tAMR\n\
             famA\tROOTX\tfamA\tCLASSA\tSUBA\tcore\tAMR\tAMR\n\
             partB\tfamB\tpartB\tCLASSB\tSUBB\tcore\tAMR\tAMR\n\
             famB\tROOTX\tfamB\tCLASSB\tSUBB\tcore\tAMR\tAMR\n\
             ROOTX\t\tROOTX\tCLASSR\tSUBR\tcore\tAMR\tAMR\n\
             relA\tfamA\trelA\tCLASSA\tSUBA\tcore\tAMR\tAMR\n",
        )
        .unwrap();
        fs::write(
            dir.join("ReferenceGeneCatalog.txt"),
            "allele\tgene_family\tproduct_name\tscope\ttype\tsubtype\tclass\tsubclass\trefseq_protein_accession\trefseq_nucleotide_accession\tgenbank_protein_accession\tgenbank_nucleotide_accession\thierarchy_node\n\
             \tpartA/partB\tfusion product\tcore\tAMR\tAMR\tCLASSA\tSUBA\tWP_F\tNG_F\t\t\tpartA,partB\n\
             relA-1\trelA\trel product\tcore\tAMR\tAMR\tCLASSA\tSUBA\tWP_R\tNG_R\t\t\trelA\n",
        )
        .unwrap();
        fs::write(
            dir.join("AMR_CDS.fa"),
            format!(
                ">WP_F|NG_F|1|2|partA|partA|fusion_product\n{part1_seq}\n\
                 >WP_F|NG_F|2|2|partB|partB|fusion_product\n{part2_seq}\n\
                 >WP_R|NG_R|1|1|relA|relA|rel_product\nTTGGCCAATTGGCCAATTGGAACCAAGGTT\n"
            ),
        )
        .unwrap();
    }

    const FUSION_SEQ: &str = "ACGTACGTAAATTTCCCGGGATATATCCCC";

    #[test]
    fn merges_identical_fusion_records_into_one_reference() {
        let dir = scratch_db_dir("fusion_merge");
        write_fusion_test_db(&dir, FUSION_SEQ, FUSION_SEQ);
        let references = load_amrfinder_references(&dir, &[ReferenceType::Amr]).unwrap();
        fs::remove_dir_all(&dir).ok();

        assert_eq!(references.len(), 2);
        let fusion = references
            .iter()
            .find(|reference| reference.protein_accession == "WP_F")
            .unwrap();
        assert_eq!(fusion.element_symbol, "partA/partB");
        assert_eq!(fusion.gene_symbol, "partA/partB");
        assert_eq!(fusion.hierarchy_node, "partA,partB");
        assert_eq!(fusion.class_name, "CLASSA");
        assert_eq!(fusion.subclass, "SUBA");
        assert_eq!(fusion.seq, FUSION_SEQ.as_bytes());
        let path_ids: Vec<&str> = fusion
            .hierarchy_path
            .iter()
            .map(|node| node.node_id.as_str())
            .collect();
        assert_eq!(path_ids, ["partA", "famA", "partB", "famB", "ROOTX"]);
    }

    #[test]
    fn keeps_fusion_parts_with_distinct_sequences() {
        let dir = scratch_db_dir("fusion_distinct");
        write_fusion_test_db(&dir, FUSION_SEQ, "GGGGGAAAAACCCCCTTTTTGGGGGAAAAA");
        let references = load_amrfinder_references(&dir, &[ReferenceType::Amr]).unwrap();
        fs::remove_dir_all(&dir).ok();

        let fusion_parts = references
            .iter()
            .filter(|reference| reference.protein_accession == "WP_F")
            .count();
        assert_eq!(fusion_parts, 2);
        assert_eq!(references.len(), 3);
    }

    #[test]
    fn merged_fusion_is_detectable() {
        use crate::detect::{DetectParams, QueryKind, detect_fasta};
        use crate::index::{IndexAlphabet, IndexBuildConfig, build_index};

        let dir = scratch_db_dir("fusion_detect");
        write_fusion_test_db(&dir, FUSION_SEQ, FUSION_SEQ);
        let references = load_amrfinder_references(&dir, &[ReferenceType::Amr]).unwrap();
        fs::remove_dir_all(&dir).ok();

        let index = build_index(
            &references,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Dna,
                k: 5,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        let fasta = format!(">contig\n{FUSION_SEQ}\n");
        let result = detect_fasta(
            &index,
            fasta.as_bytes(),
            "sample",
            QueryKind::Direct,
            &DetectParams {
                min_gene_fraction: 0.5,
                ..DetectParams::default()
            },
        )
        .unwrap();
        assert!(
            result
                .hits
                .iter()
                .any(|hit| hit.element_symbol.as_deref() == Some("partA/partB")
                    && hit.call_type == "gene"),
            "fusion gene was not reported: {:?}",
            result.hits
        );
    }

    #[test]
    fn merge_catalog_entry_rejects_biological_conflicts() {
        let mut existing = CatalogEntry::from_row(&sample_catalog_row("WP_000649751.1")).unwrap();
        let incoming = CatalogEntry::from_row(&CatalogRow {
            hierarchy_node: "different_node",
            ..sample_catalog_row("WP_000649751.1")
        })
        .unwrap();
        let err = merge_catalog_entry(&mut existing, &incoming).unwrap_err();
        assert!(err.to_string().contains("existing="));
    }
}
