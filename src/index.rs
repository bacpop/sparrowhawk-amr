use crate::amrfinder_db::{AmrReference, HierarchyNode};
use crate::kmer::{DnaKmerIter, ProteinKmerIter};
use anyhow::{Context, ensure};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
#[cfg(not(target_family = "wasm"))]
use std::fs;
#[cfg(not(target_family = "wasm"))]
use std::path::Path;

const INDEX_MAGIC: &[u8; 8] = b"SHAMR004";

pub type GeneId = u16;
pub type UnitId = u16;
pub type StringId = u16;

/// This is the string subindex
#[derive(Debug, Default)]
struct StringInterner {
    strings: Vec<String>,
    map: HashMap<String, StringId>,
}

impl StringInterner {
    fn intern(&mut self, value: &str) -> anyhow::Result<StringId> {
        // See if it's inside...
        if let Some(&id) = self.map.get(value) {
            return Ok(id);
        }

        // if not, add it if possible
        let id = StringId::try_from(self.strings.len())
            .context("too many strings for u16 string IDs")?;
        self.strings.push(value.to_string());
        self.map.insert(value.to_string(), id);
        Ok(id)
    }
}

fn checked_gene_id(value: usize) -> anyhow::Result<GeneId> {
    GeneId::try_from(value).context("too many genes for u16 gene IDs")
}

fn checked_unit_id(value: usize) -> anyhow::Result<UnitId> {
    UnitId::try_from(value).context("too many report units for u16 unit IDs")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReportUnitKind {
    ExactGene,
    HierarchyNode,
}

impl ReportUnitKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExactGene => "exact_gene",
            Self::HierarchyNode => "hierarchy_node",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum IndexAlphabet {
    Dna,
    Protein,
}

impl IndexAlphabet {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Dna => "dna",
            Self::Protein => "protein",
        }
    }
}

impl Default for IndexAlphabet {
    fn default() -> Self {
        Self::Dna
    }
}


// Struct to work with per-gene info: this is temporal, but used for debugging!
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneEntry {
    pub id: StringId,
    pub element_symbol: StringId,
    pub gene_symbol: StringId,
    pub allele_symbol: StringId,
    pub protein_accession: StringId,
    pub nucleotide_accession: StringId,
    pub gene_group: StringId,
    pub class_name: StringId,
    pub subclass: StringId,
    pub type_name: StringId,
    pub subtype: StringId,
    pub hierarchy_node: StringId,
    pub product: StringId,
    pub length: usize,
    pub gene_specific_kmers: usize,
    pub report_unit_id: UnitId,
    pub exact_unit_eligible: bool,
}


// Struct that will become part of the index
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportUnit {
    pub id: StringId,
    pub label: StringId,
    pub gene_id: Option<GeneId>,
    pub element_symbol: Option<StringId>,
    pub gene_symbol: Option<StringId>,
    pub allele_symbol: Option<StringId>,
    pub gene_group: StringId,
    pub hierarchy_node: StringId,
    pub class_name: StringId,
    pub subclass: StringId,
    pub type_name: StringId,
    pub subtype: StringId,
    pub product: StringId,
    pub member_count: usize,
    pub member_gene_ids: Vec<GeneId>,
    pub diagnostic_kmers: usize,
    pub weak: bool,
    pub ancestor_unit_ids: Vec<UnitId>,
}

// The index itself
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmrIndex {
    pub magic: [u8; 8],
    pub db_version: String,
    pub alphabet: IndexAlphabet,
    pub k: usize,
    pub min_exact_gene_kmers: usize,
    pub min_hierarchy_unit_kmers: usize,
    pub strings: Vec<String>,
    pub genes: Vec<GeneEntry>,
    pub units: Vec<ReportUnit>,
    pub kmer_codes: Vec<u64>,
    pub unit_ids: Vec<UnitId>,
}

#[derive(Debug, Clone)]
pub struct IndexBuildConfig {
    pub alphabet: IndexAlphabet,
    pub k: usize,
    pub min_exact_gene_kmers: usize,
    pub min_hierarchy_unit_kmers: usize,
}

impl Default for IndexBuildConfig {
    fn default() -> Self {
        Self {
            alphabet: IndexAlphabet::Dna,
            k: 31,
            min_exact_gene_kmers: 20,
            min_hierarchy_unit_kmers: 20,
        }
    }
}

impl ReportUnit {
    pub fn kind(&self) -> ReportUnitKind {
        if self.gene_id.is_some() {
            ReportUnitKind::ExactGene
        } else {
            ReportUnitKind::HierarchyNode
        }
    }

    pub fn call_type(&self) -> &'static str {
        match self.kind() {
            ReportUnitKind::ExactGene => "gene",
            ReportUnitKind::HierarchyNode => "gene_group",
        }
    }
}

impl AmrIndex {
    pub fn lookup(&self, code: u64) -> Option<usize> {
        let idx = self.kmer_codes.binary_search(&code).ok()?;
        Some(self.unit_ids[idx] as usize)
    }

    pub fn unit_specific_kmers(&self, unit_id: usize) -> HashSet<u64> {
        self.kmer_codes
            .iter()
            .zip(&self.unit_ids)
            .filter_map(|(&kmer, &assigned)| (assigned as usize == unit_id).then_some(kmer))
            .collect()
    }

    pub fn string(&self, id: StringId) -> &str {
        &self.strings[id as usize]
    }

    pub fn optional_string(&self, id: Option<StringId>) -> Option<String> {
        id.map(|value| self.string(value).to_string())
            .filter(|value| !value.is_empty())
    }

    /// Debugging/info method
    pub fn stats_string(&self) -> String {
        let exact_units = self
            .units
            .iter()
            .filter(|unit| unit.kind() == ReportUnitKind::ExactGene)
            .count();
        let hierarchy_units = self.units.len().saturating_sub(exact_units);
        let collapsed_genes = self
            .genes
            .iter()
            .filter(|gene| !gene.exact_unit_eligible)
            .count();
        let genes_at_or_below_exact_threshold = self
            .genes
            .iter()
            .filter(|gene| gene.gene_specific_kmers <= self.min_exact_gene_kmers)
            .count();
        let weak_hierarchy_units = self
            .units
            .iter()
            .filter(|unit| unit.kind() == ReportUnitKind::HierarchyNode && unit.weak)
            .count();
        let type_counts = format_counts(self.genes.iter().map(|gene| self.string(gene.type_name)));
        let subtype_counts = format_counts(self.genes.iter().map(|gene| self.string(gene.subtype)));
        format!(
            "db_version={}\nalphabet={}\nk={}\nmin_exact_gene_kmers={}\nmin_hierarchy_unit_kmers={}\ngenes={}\ngenes_at_or_below_exact_threshold={}\ncollapsed_genes={}\nreport_units={}\nexact_gene_units={}\nhierarchy_units={}\nweak_hierarchy_units={}\nretained_kmers={}\nunit_specific_kmers={}\nindex_bytes_uncompressed_estimate={}\ntype_counts={}\nsubtype_counts={}\n",
            self.db_version,
            self.alphabet.as_str(),
            self.k,
            self.min_exact_gene_kmers,
            self.min_hierarchy_unit_kmers,
            self.genes.len(),
            genes_at_or_below_exact_threshold,
            collapsed_genes,
            self.units.len(),
            exact_units,
            hierarchy_units,
            weak_hierarchy_units,
            self.kmer_codes.len(),
            self.units
                .iter()
                .map(|unit| unit.diagnostic_kmers)
                .sum::<usize>(),
            self.estimated_bytes(),
            type_counts,
            subtype_counts,
        )
    }

    fn estimated_bytes(&self) -> usize {
        self.kmer_codes.len() * 8 + self.unit_ids.len() * 2
    }
}


fn format_counts<'a>(values: impl Iterator<Item = &'a str>) -> String {
    let mut counts = BTreeMap::<String, usize>::new();
    for value in values {
        if !value.is_empty() {
            *counts.entry(value.to_string()).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .map(|(key, count)| format!("{key}:{count}"))
        .collect::<Vec<_>>()
        .join(",")
}


/// Constructs the index
pub fn build_index(
    references: &[AmrReference],
    config: &IndexBuildConfig,
) -> anyhow::Result<AmrIndex> {

    // First, initsssss
    match config.alphabet {
        IndexAlphabet::Dna => ensure!(
            (1..=31).contains(&config.k),
            "DNA k must be between 1 and 31"
        ),
        IndexAlphabet::Protein => ensure!(
            (1..=12).contains(&config.k),
            "protein k must be between 1 and 12"
        ),
    }

    let mut interner = StringInterner::default();
    let mut genes = Vec::<GeneEntry>::new();
    let mut gene_paths = Vec::<Vec<String>>::new();
    let mut node_meta = HashMap::<String, NodeMeta>::new();
    let mut node_to_genes = HashMap::<String, HashSet<usize>>::new();
    let mut raw_kmers = HashMap::<u64, Vec<usize>>::new();

    // Parse all info from references into geneentries for the debug info of the index, and get their kmers
    for reference in references {
        let gene_id = genes.len();
        checked_gene_id(gene_id)?; // does it fit in a u16?
        let gene_label = if reference.allele_symbol.is_empty() {
            reference.gene_symbol.clone()
        } else {
            reference.allele_symbol.clone()
        };
        let path = normalized_path(reference); // Just create one if there is none
        for node in &path {
            node_meta
                .entry(node.node_id.clone())
                .or_insert_with(|| NodeMeta::from_node(node));
            node_to_genes
                .entry(node.node_id.clone())
                .or_default()
                .insert(gene_id);
        }
        let path_ids: Vec<String> = path.into_iter().map(|node| node.node_id).collect();
        genes.push(GeneEntry {
            id: interner.intern(&format!("{}|{}", gene_label, reference.protein_accession))?,
            element_symbol: interner.intern(&reference.element_symbol)?,
            gene_symbol: interner.intern(&reference.gene_symbol)?,
            allele_symbol: interner.intern(&reference.allele_symbol)?,
            protein_accession: interner.intern(&reference.protein_accession)?,
            nucleotide_accession: interner.intern(&reference.nucleotide_accession)?,
            gene_group: interner.intern(&reference.family)?,
            class_name: interner.intern(&reference.class_name)?,
            subclass: interner.intern(&reference.subclass)?,
            type_name: interner.intern(&reference.type_name)?,
            subtype: interner.intern(&reference.subtype)?,
            hierarchy_node: interner.intern(&reference.hierarchy_node)?,
            product: interner.intern(&reference.product)?,
            length: reference.seq.len(),
            gene_specific_kmers: 0,
            report_unit_id: 0,
            exact_unit_eligible: false,
        });
        gene_paths.push(path_ids);

        // This could be done better...
        let Some(iter) = reference_kmers(&reference.seq, config.k, config.alphabet) else {
            continue;
        };
        let mut seen_in_gene = HashSet::new();
        for kmer in iter {
            if seen_in_gene.insert(kmer) {
                raw_kmers.entry(kmer).or_default().push(gene_id);
            }
        }
    }


    // Now, let's start getting values for setting the report units.
    // First, gene unique counts
    let mut gene_unique_counts = vec![0usize; genes.len()];
    for refs in raw_kmers.values() {
        if refs.len() == 1 {
            gene_unique_counts[refs[0]] += 1;
        }
    }

    // See which genes have enough unique k-mers
    let exact_eligible: Vec<bool> = gene_unique_counts
        .iter()
        .map(|&count| count > config.min_exact_gene_kmers)
        .collect();
    for (gene, (&count, &eligible)) in genes
        .iter_mut()
        .zip(gene_unique_counts.iter().zip(&exact_eligible))
    {
        gene.gene_specific_kmers = count;
        gene.exact_unit_eligible = eligible;
    }

    // Now, let's go for the report units that cannot be unique genes/alleles
    // This is a map that gets you the counts per node
    let node_candidate_counts =
        hierarchy_candidate_counts(&raw_kmers, &gene_paths, &exact_eligible);

    // selection magic
    let selected_node_by_gene = select_hierarchy_units(
        &genes,
        &gene_paths,
        &exact_eligible,
        &node_candidate_counts,
        config.min_hierarchy_unit_kmers,
    );

    // With this, we can construct the units
    let mut units = Vec::<ReportUnit>::new();
    let mut node_unit_ids = HashMap::<String, UnitId>::new();
    // First, exact genes
    for gene_id in 0..genes.len() {
        if !exact_eligible[gene_id] {
            continue;
        }
        let unit_id = checked_unit_id(units.len())?;
        genes[gene_id].report_unit_id = unit_id;
        units.push(exact_unit(checked_gene_id(gene_id)?, &genes[gene_id]));
    }

    // Now the superior hiearchical nodes
    let selected_nodes: BTreeSet<String> = selected_node_by_gene
        .iter()
        .filter_map(|node| node.clone())
        .collect();
    for node_id in selected_nodes {
        let unit_id = checked_unit_id(units.len())?;
        node_unit_ids.insert(node_id.clone(), unit_id);
        units.push(hierarchy_unit(
            &node_id,
            &node_meta,
            &node_to_genes,
            &node_candidate_counts,
            config.min_hierarchy_unit_kmers,
            &mut interner,
        )?);
    }

    // Assign an exact gene (for debug) to a report unit
    for (gene_id, selected_node) in selected_node_by_gene.iter().enumerate() {
        if exact_eligible[gene_id] {
            continue;
        }
        let Some(node_id) = selected_node else {
            continue;
        };
        let Some(&unit_id) = node_unit_ids.get(node_id) else {
            continue;
        };
        genes[gene_id].report_unit_id = unit_id;
    }

    // fill all ancestors, this is used to not report them once a more specific call has been made
    fill_unit_ancestors(&mut units, &gene_paths, &node_unit_ids, &interner.strings);

    // Write the k-mers
    let mut retained = Vec::<(u64, UnitId)>::new();
    for (kmer, refs) in raw_kmers {
        if refs.len() == 1 && exact_eligible[refs[0]] {
            retained.push((kmer, genes[refs[0]].report_unit_id));
            continue;
        }
        let Some(unit_id) = lowest_selected_hierarchy_unit(&refs, &gene_paths, &node_unit_ids)
        else {
            continue;
        };
        retained.push((kmer, unit_id));
    }
    retained.sort_by_key(|(kmer, _)| *kmer);

    for &(_, unit_id) in &retained {
        if let Some(unit) = units.get_mut(unit_id as usize) {
            unit.diagnostic_kmers += 1;
        }
    }

    // Finally, the db version
    let db_version = references
        .first()
        .map(|reference| reference.db_version.clone())
        .unwrap_or_else(|| "unknown".to_string());

    Ok(AmrIndex {
        magic: *INDEX_MAGIC,
        db_version,
        alphabet: config.alphabet,
        k: config.k,
        min_exact_gene_kmers: config.min_exact_gene_kmers,
        min_hierarchy_unit_kmers: config.min_hierarchy_unit_kmers,
        strings: interner.strings,
        genes,
        units,
        kmer_codes: retained.iter().map(|(kmer, _)| *kmer).collect(),
        unit_ids: retained.iter().map(|(_, unit_id)| *unit_id).collect(),
    })
}


#[cfg(not(target_family = "wasm"))]
pub fn save_index(index: &AmrIndex, path: &Path) -> anyhow::Result<()> {
    let bytes = bincode::serialize(index).context("serialize AMR index")?;
    fs::write(path, bytes).with_context(|| format!("write index {}", path.display()))
}

#[cfg(not(target_family = "wasm"))]
pub fn load_index(path: &Path) -> anyhow::Result<AmrIndex> {
    let bytes = fs::read(path).with_context(|| format!("read index {}", path.display()))?;
    load_index_from_bytes(&bytes)
}

pub fn load_index_from_bytes(bytes: &[u8]) -> anyhow::Result<AmrIndex> {
    let index: AmrIndex = bincode::deserialize(&bytes).context("deserialize AMR index")?;
    ensure!(index.magic == *INDEX_MAGIC, "unsupported AMR index format");
    match index.alphabet {
        IndexAlphabet::Dna => ensure!((1..=31).contains(&index.k), "invalid DNA k in AMR index"),
        IndexAlphabet::Protein => ensure!(
            (1..=12).contains(&index.k),
            "invalid protein k in AMR index"
        ),
    }
    ensure!(
        index.kmer_codes.len() == index.unit_ids.len(),
        "corrupt AMR index assignment arrays"
    );
    ensure!(
        index
            .unit_ids
            .iter()
            .all(|&unit_id| (unit_id as usize) < index.units.len()),
        "corrupt AMR index unit ids"
    );
    Ok(index)
}

fn reference_kmers(seq: &[u8], k: usize, alphabet: IndexAlphabet) -> Option<Vec<u64>> {
    match alphabet {
        IndexAlphabet::Dna => {
            DnaKmerIter::new(seq, k).map(|iter| iter.map(|(_pos, kmer)| kmer).collect())
        }
        IndexAlphabet::Protein => {
            ProteinKmerIter::new(seq, k).map(|iter| iter.map(|(_pos, kmer)| kmer).collect())
        }
    }
}


#[derive(Debug, Clone)]
struct NodeMeta {
    label: String,
    class_name: String,
    subclass: String,
    type_name: String,
    subtype: String,
}

impl NodeMeta {
    fn from_node(node: &HierarchyNode) -> Self {
        Self {
            label: if node.symbol.is_empty() {
                node.node_id.clone()
            } else {
                node.symbol.clone()
            },
            class_name: node.class_name.clone(),
            subclass: node.subclass.clone(),
            type_name: node.type_name.clone(),
            subtype: node.subtype.clone(),
        }
    }
}

fn normalized_path(reference: &AmrReference) -> Vec<HierarchyNode> {
    if !reference.hierarchy_path.is_empty() {
        return reference.hierarchy_path.clone();
    }
    let node_id = first_non_empty([
        reference.hierarchy_node.as_str(),
        reference.family.as_str(),
        reference.element_symbol.as_str(),
    ]);
    vec![HierarchyNode {
        node_id: node_id.to_string(),
        parent_node_id: String::new(),
        symbol: node_id.to_string(),
        class_name: reference.class_name.clone(),
        subclass: reference.subclass.clone(),
        scope: reference.scope.clone(),
        type_name: reference.type_name.clone(),
        subtype: reference.subtype.clone(),
        reportable: reference.reportable,
    }]
}


/// Get the counts for all nodes that might be susceptible of giving a report unit over a single gene/allele
/// (i.e. ignoring those that we can call directly).
fn hierarchy_candidate_counts(
    raw_kmers: &HashMap<u64, Vec<usize>>,
    gene_paths: &[Vec<String>],
    exact_eligible: &[bool],
) -> HashMap<String, usize> {
    let mut counts = HashMap::<String, usize>::new();
    for refs in raw_kmers.values() {
        if refs.len() == 1 && exact_eligible[refs[0]] {
            continue;
        }
        for node_id in common_path_nodes(refs, gene_paths) {
            *counts.entry(node_id).or_default() += 1;
        }
    }
    counts
}


/// Main function for selecting which other superior hierarchy units we choose
fn select_hierarchy_units(
    genes: &[GeneEntry],
    gene_paths: &[Vec<String>],
    exact_eligible: &[bool],
    node_candidate_counts: &HashMap<String, usize>,
    min_hierarchy_unit_kmers: usize,
) -> Vec<Option<String>> {
    let mut selected = vec![None; genes.len()];
    for gene_id in 0..genes.len() {
        if exact_eligible[gene_id] {
            continue;
        }

        let path = &gene_paths[gene_id];

        selected[gene_id] = path
            .iter()
            .find(|node_id| {
                node_candidate_counts
                    .get(*node_id)
                    .copied()
                    .unwrap_or_default()
                    >= min_hierarchy_unit_kmers
            })
            .cloned()
            .or_else(|| { // If not a suitable one has been found, let's get the largest one possible
                path.iter()
                    .max_by_key(|node_id| {
                        node_candidate_counts
                            .get(*node_id)
                            .copied()
                            .unwrap_or_default()
                    })
                    .cloned()
            });
    }

    selected
}

fn exact_unit(gene_id: GeneId, gene: &GeneEntry) -> ReportUnit {
    ReportUnit {
        id: gene.id,
        label: gene.element_symbol,
        gene_id: Some(gene_id),
        element_symbol: Some(gene.element_symbol),
        gene_symbol: Some(gene.gene_symbol),
        allele_symbol: Some(gene.allele_symbol),
        gene_group: gene.gene_group,
        hierarchy_node: gene.hierarchy_node,
        class_name: gene.class_name,
        subclass: gene.subclass,
        type_name: gene.type_name,
        subtype: gene.subtype,
        product: gene.product,
        member_count: 1,
        member_gene_ids: vec![gene_id],
        diagnostic_kmers: 0,
        weak: false,
        ancestor_unit_ids: Vec::new(),
    }
}

fn hierarchy_unit(
    node_id: &str,
    node_meta: &HashMap<String, NodeMeta>,
    node_to_genes: &HashMap<String, HashSet<usize>>,
    node_candidate_counts: &HashMap<String, usize>,
    min_hierarchy_unit_kmers: usize,
    interner: &mut StringInterner,
) -> anyhow::Result<ReportUnit> {
    let meta = node_meta.get(node_id);
    let mut member_gene_ids: Vec<GeneId> = node_to_genes
        .get(node_id)
        .into_iter()
        .flat_map(|genes| genes.iter())
        .map(|&gene_id| checked_gene_id(gene_id))
        .collect::<anyhow::Result<Vec<_>>>()?;
    member_gene_ids.sort_unstable();
    let candidate_count = node_candidate_counts
        .get(node_id)
        .copied()
        .unwrap_or_default();
    Ok(ReportUnit {
        id: interner.intern(node_id)?,
        label: interner.intern(meta.map(|meta| meta.label.as_str()).unwrap_or(node_id))?,
        gene_id: None,
        element_symbol: None,
        gene_symbol: None,
        allele_symbol: None,
        gene_group: interner.intern(node_id)?,
        hierarchy_node: interner.intern(node_id)?,
        class_name: interner.intern(
            meta.map(|meta| meta.class_name.as_str())
                .unwrap_or_default(),
        )?,
        subclass: interner.intern(meta.map(|meta| meta.subclass.as_str()).unwrap_or_default())?,
        type_name: interner.intern(meta.map(|meta| meta.type_name.as_str()).unwrap_or_default())?,
        subtype: interner.intern(meta.map(|meta| meta.subtype.as_str()).unwrap_or_default())?,
        product: interner.intern("")?,
        member_count: member_gene_ids.len(),
        member_gene_ids,
        diagnostic_kmers: 0,
        weak: candidate_count < min_hierarchy_unit_kmers,
        ancestor_unit_ids: Vec::new(),
    })
}

// This helper collects all ancestors of a unit report
fn fill_unit_ancestors(
    units: &mut [ReportUnit],
    gene_paths: &[Vec<String>],
    node_unit_ids: &HashMap<String, UnitId>,
    strings: &[String],
) {
    for unit_id in 0..units.len() {
        let ancestor_nodes: Vec<&String> = match units[unit_id].kind() {
            ReportUnitKind::ExactGene => {
                let Some(gene_id) = units[unit_id].gene_id else {
                    continue;
                };
                gene_paths[gene_id as usize].iter().collect()
            }
            ReportUnitKind::HierarchyNode => {
                let Some(gene_id) = units[unit_id].member_gene_ids.first() else {
                    continue;
                };
                let Some(position) = gene_paths[*gene_id as usize].iter().position(|node_id| {
                    node_id == &strings[units[unit_id].hierarchy_node as usize]
                }) else {
                    continue;
                };
                gene_paths[*gene_id as usize][position + 1..]
                    .iter()
                    .collect()
            }
        };
        let mut ancestors = Vec::new();
        for node_id in ancestor_nodes {
            let Some(&ancestor_unit_id) = node_unit_ids.get(node_id) else {
                continue;
            };
            if ancestor_unit_id as usize != unit_id && !ancestors.contains(&ancestor_unit_id) {
                ancestors.push(ancestor_unit_id);
            }
        }
        units[unit_id].ancestor_unit_ids = ancestors;
    }
}


fn lowest_selected_hierarchy_unit(
    refs: &[usize],
    gene_paths: &[Vec<String>],
    node_unit_ids: &HashMap<String, UnitId>,
) -> Option<UnitId> {
    common_path_nodes(refs, gene_paths)
        .into_iter()
        .find_map(|node_id| node_unit_ids.get(&node_id).copied())
}


/// S
fn common_path_nodes(refs: &[usize], gene_paths: &[Vec<String>]) -> Vec<String> {
    let Some((&first, rest)) = refs.split_first() else {
        return Vec::new();
    }; // Just a check for empty slices

    gene_paths[first]
        .iter()
        .filter(|node_id| {
            rest.iter()
                .all(|&gene_id| gene_paths[gene_id].iter().any(|other| other == *node_id))
        })
        .cloned()
        .collect()
}

fn first_non_empty<'a>(values: [&'a str; 3]) -> &'a str {
    values
        .into_iter()
        .find(|value| !value.is_empty())
        .unwrap_or("")
}



// =============================================== TESTs

#[cfg(test)]
mod tests {
    use super::*;

    fn reference(id: &str, node: &str, seq: &[u8]) -> AmrReference {
        AmrReference {
            protein_accession: id.to_string(),
            nucleotide_accession: String::new(),
            element_symbol: id.to_string(),
            gene_symbol: id.to_string(),
            allele_symbol: id.to_string(),
            product: String::new(),
            family: node.to_string(),
            class_name: String::new(),
            subclass: String::new(),
            hierarchy_node: String::new(),
            scope: "core".to_string(),
            type_name: "AMR".to_string(),
            subtype: "AMR".to_string(),
            reportable: 2,
            hierarchy_path: vec![HierarchyNode {
                node_id: node.to_string(),
                parent_node_id: String::new(),
                symbol: node.to_string(),
                class_name: String::new(),
                subclass: String::new(),
                scope: "core".to_string(),
                type_name: "AMR".to_string(),
                subtype: "AMR".to_string(),
                reportable: 2,
            }],
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
        assert_eq!(index.genes.len(), 2);
        assert_eq!(index.units.len(), 2);
        assert!(!index.kmer_codes.is_empty());
        assert_eq!(index.kmer_codes.len(), index.unit_ids.len());
    }

    #[test]
    fn exact_threshold_is_strictly_greater_than_configured_value() {
        let refs = vec![reference("g1", "node", b"ACGTACGTAC")];
        let loose = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Dna,
                k: 5,
                min_exact_gene_kmers: 0,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        assert!(loose.genes[0].exact_unit_eligible);

        let exact_count = loose.genes[0].gene_specific_kmers;
        let strict = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Dna,
                k: 5,
                min_exact_gene_kmers: exact_count,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        assert!(!strict.genes[0].exact_unit_eligible);
        assert_eq!(
            strict.units[strict.genes[0].report_unit_id as usize].kind(),
            ReportUnitKind::HierarchyNode
        );
    }

    #[test]
    fn protein_index_uses_hierarchy_collapse_for_weak_exact_genes() {
        let refs = vec![
            reference("p1", "node", b"MAAAAAAAK"),
            reference("p2", "node", b"MAAAAAAAR"),
        ];
        let index = build_index(
            &refs,
            &IndexBuildConfig {
                alphabet: IndexAlphabet::Protein,
                k: 3,
                min_exact_gene_kmers: 5,
                min_hierarchy_unit_kmers: 1,
            },
        )
        .unwrap();
        assert_eq!(index.alphabet, IndexAlphabet::Protein);
        assert!(index.genes.iter().all(|gene| !gene.exact_unit_eligible));
        assert_eq!(index.units.len(), 1);
        assert_eq!(index.units[0].kind(), ReportUnitKind::HierarchyNode);
        assert!(!index.kmer_codes.is_empty());
    }
    #[test]
    fn index_preserves_reference_category_metadata() {
        let mut refs = vec![reference("stress1", "metal_node", b"ACGTACGTAC")];
        refs[0].type_name = "STRESS".to_string();
        refs[0].subtype = "METAL".to_string();
        refs[0].hierarchy_path[0].type_name = "STRESS".to_string();
        refs[0].hierarchy_path[0].subtype = "METAL".to_string();
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
        assert_eq!(index.string(index.genes[0].type_name), "STRESS");
        assert_eq!(index.string(index.genes[0].subtype), "METAL");
        let unit = &index.units[index.genes[0].report_unit_id as usize];
        assert_eq!(index.string(unit.type_name), "STRESS");
        assert_eq!(index.string(unit.subtype), "METAL");
        assert!(index.stats_string().contains("type_counts=STRESS:1"));
    }
}
