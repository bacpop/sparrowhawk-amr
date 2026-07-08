#[cfg(target_family = "wasm")]
fn main() {}

#[cfg(not(target_family = "wasm"))]
mod native {
    use anyhow::{Context, bail};
    use clap::{Parser, Subcommand, ValueEnum};
    use sparrowhawk_amr::{
        DebugMissesConfig, DetectParams, GeneCallerConfig, IndexAlphabet, IndexBuildConfig,
        QueryKind, ReferenceType, RefinementMode, TruthKmerEvidenceConfig, build_index,
        debug_amrfinder_misses, detect_fasta, detect_protein_fasta,
        load_amrfinder_protein_references, load_amrfinder_references, load_index, run_gene_caller,
        save_index,
    };
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::process::Command;

    const DEFAULT_DB_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest";
    const DB_FILES: &[&str] = &[
        "AMR_CDS.fa",
        "AMRProt.fa",
        "ReferenceGeneCatalog.txt",
        "ReferenceGeneHierarchy.txt",
        "fam.tsv",
        "version.txt",
    ];

    #[derive(Parser)]
    #[command(author, version, about = "Offline AMRFinderPlus k-mer AMR detector")]
    struct Cli {
        #[command(subcommand)]
        command: CommandKind,
    }

    #[derive(Subcommand)]
    enum CommandKind {
        Db {
            #[command(subcommand)]
            command: DbCommand,
        },
        Index {
            #[command(subcommand)]
            command: IndexCommand,
        },
        Detect {
            #[command(subcommand)]
            command: DetectCommand,
        },
        Genes {
            #[command(subcommand)]
            command: GenesCommand,
        },
        Eval {
            #[command(subcommand)]
            command: EvalCommand,
        },
    }

    #[derive(Subcommand)]
    enum DbCommand {
        Fetch {
            #[arg(long)]
            out_dir: PathBuf,
            #[arg(long, default_value = DEFAULT_DB_URL)]
            base_url: String,
        },
    }

    #[derive(Subcommand)]
    enum IndexCommand {
        Build {
            #[arg(long)]
            db_dir: PathBuf,
            #[arg(long)]
            out: PathBuf,
            #[arg(long, value_enum, default_value_t = IndexAlphabetArg::Dna)]
            alphabet: IndexAlphabetArg,
            #[arg(long)]
            k: Option<usize>,
            #[arg(long)]
            min_exact_gene_kmers: Option<usize>,
            #[arg(long)]
            min_hierarchy_unit_kmers: Option<usize>,
            #[arg(long, value_enum, value_delimiter = ',', default_values_t = [
                ReferenceTypeArg::Amr,
                ReferenceTypeArg::Stress,
                ReferenceTypeArg::Virulence,
            ])]
            include_types: Vec<ReferenceTypeArg>,
        },
        Stats {
            #[arg(long)]
            index: PathBuf,
        },
        ReportMap {
            #[arg(long)]
            index: PathBuf,
            #[arg(long)]
            out: Option<PathBuf>,
        },
        UnitStats {
            #[arg(long)]
            index: PathBuf,
            #[arg(long)]
            db_dir: Option<PathBuf>,
            #[arg(long)]
            out: Option<PathBuf>,
        },
    }

    #[derive(Subcommand)]
    enum DetectCommand {
        Direct {
            #[arg(long)]
            index: PathBuf,
            #[arg(long)]
            fasta: PathBuf,
            #[arg(long)]
            sample_name: Option<String>,
            #[command(flatten)]
            params: DetectArgs,
        },
        Cds {
            #[arg(long)]
            index: PathBuf,
            #[arg(long)]
            assembly: Option<PathBuf>,
            #[arg(long)]
            cds_fasta: Option<PathBuf>,
            #[arg(long)]
            protein: bool,
            #[arg(long)]
            protein_fasta: Option<PathBuf>,
            #[arg(long, default_value = "gene_calls")]
            out_dir: PathBuf,
            #[command(flatten)]
            orphos: OrphosArgs,
            #[arg(long)]
            sample_name: Option<String>,
            #[command(flatten)]
            params: DetectArgs,
        },
    }

    #[derive(Subcommand)]
    enum GenesCommand {
        Call {
            #[arg(long)]
            assembly: PathBuf,
            #[arg(long)]
            out_dir: PathBuf,
            #[command(flatten)]
            orphos: OrphosArgs,
            #[arg(long)]
            sample_name: Option<String>,
        },
    }

    #[derive(Subcommand)]
    enum EvalCommand {
        NativeAmrfinder {
            #[arg(long)]
            amrfinder_path: Option<PathBuf>,
            #[arg(long)]
            db_dir: PathBuf,
            #[arg(long)]
            assembly: PathBuf,
            #[arg(long)]
            out: PathBuf,
        },
        DebugAmrfinderMisses {
            #[arg(long)]
            index: PathBuf,
            #[arg(long)]
            assembly: PathBuf,
            #[arg(long)]
            amrfinder_tsv: PathBuf,
            #[arg(long)]
            detector_json: PathBuf,
            #[arg(long)]
            db_dir: Option<PathBuf>,
            #[arg(long, default_value_t = 17)]
            refinement_k: usize,
            #[arg(long, default_value_t = 5)]
            missing_kmer_limit: usize,
            #[arg(long)]
            out: Option<PathBuf>,
        },
        TruthKmerEvidence {
            #[arg(long)]
            index: PathBuf,
            #[arg(long)]
            assembly: PathBuf,
            #[arg(long)]
            amrfinder_tsv: PathBuf,
            #[arg(long)]
            detector_json: PathBuf,
            #[arg(long)]
            db_dir: Option<PathBuf>,
            #[arg(long, value_enum, value_delimiter = ',', default_values_t = [
                ReferenceTypeArg::Amr,
                ReferenceTypeArg::Stress,
                ReferenceTypeArg::Virulence,
            ])]
            include_types: Vec<ReferenceTypeArg>,
            #[arg(long, default_value_t = 0.10)]
            min_gene_fraction: f64,
            #[arg(long, default_value_t = 0.10)]
            min_family_fraction: f64,
            #[arg(long)]
            out: Option<PathBuf>,
        },
    }

    #[derive(Parser, Debug, Clone)]
    struct DetectArgs {
        #[arg(long, default_value_t = 0.10)]
        min_gene_fraction: f64,
        #[arg(
            long = "min-gene-group-fraction",
            alias = "min-family-fraction",
            default_value_t = 0.10
        )]
        min_gene_group_fraction: f64,
        #[arg(long, default_value_t = 0.01)]
        seed_gene_fraction: f64,
        #[arg(long, default_value_t = 3)]
        seed_gene_hits: usize,
        #[arg(long, value_enum, default_value_t = RefinementArg::None)]
        refinement_mode: RefinementArg,
        #[arg(long, default_value_t = 21)]
        refinement_k: usize,
    }

    #[derive(Debug, Clone, Copy, ValueEnum)]
    enum RefinementArg {
        None,
        Split,
        Lowk,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
    enum IndexAlphabetArg {
        Dna,
        Protein,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
    enum ReferenceTypeArg {
        Amr,
        Stress,
        Virulence,
    }

    #[derive(Parser, Debug, Clone, Default)]
    struct OrphosArgs {
        #[arg(long)]
        orphos_metagenomic: bool,
        #[arg(long)]
        orphos_closed_ends: bool,
        #[arg(long)]
        orphos_mask_n_runs: bool,
        #[arg(long)]
        orphos_force_non_sd: bool,
        #[arg(long)]
        orphos_translation_table: Option<u8>,
    }

    impl From<RefinementArg> for RefinementMode {
        fn from(value: RefinementArg) -> Self {
            match value {
                RefinementArg::None => Self::None,
                RefinementArg::Split => Self::Split,
                RefinementArg::Lowk => Self::LowK,
            }
        }
    }

    impl From<IndexAlphabetArg> for IndexAlphabet {
        fn from(value: IndexAlphabetArg) -> Self {
            match value {
                IndexAlphabetArg::Dna => Self::Dna,
                IndexAlphabetArg::Protein => Self::Protein,
            }
        }
    }

    impl From<ReferenceTypeArg> for ReferenceType {
        fn from(value: ReferenceTypeArg) -> Self {
            match value {
                ReferenceTypeArg::Amr => Self::Amr,
                ReferenceTypeArg::Stress => Self::Stress,
                ReferenceTypeArg::Virulence => Self::Virulence,
            }
        }
    }

    impl From<DetectArgs> for DetectParams {
        fn from(value: DetectArgs) -> Self {
            Self {
                min_gene_fraction: value.min_gene_fraction,
                min_gene_group_fraction: value.min_gene_group_fraction,
                seed_gene_fraction: value.seed_gene_fraction,
                seed_gene_hits: value.seed_gene_hits,
                refinement_mode: value.refinement_mode.into(),
                refinement_k: value.refinement_k,
            }
        }
    }

    pub fn run() -> anyhow::Result<()> {
        let cli = Cli::parse();
        match cli.command {
            CommandKind::Db { command } => match command {
                DbCommand::Fetch { out_dir, base_url } => fetch_db(&out_dir, &base_url),
            },
            CommandKind::Index { command } => match command {
                IndexCommand::Build {
                    db_dir,
                    out,
                    alphabet,
                    k,
                    min_exact_gene_kmers,
                    min_hierarchy_unit_kmers,
                    include_types,
                } => {
                    let alphabet: IndexAlphabet = alphabet.into();
                    let k = k.unwrap_or(match alphabet {
                        IndexAlphabet::Dna => 31,
                        IndexAlphabet::Protein => 5,
                    });
                    let min_exact_gene_kmers = min_exact_gene_kmers.unwrap_or(match alphabet {
                        IndexAlphabet::Dna => 20,
                        IndexAlphabet::Protein => 5,
                    });
                    let min_hierarchy_unit_kmers =
                        min_hierarchy_unit_kmers.unwrap_or(match alphabet {
                            IndexAlphabet::Dna => 20,
                            IndexAlphabet::Protein => 5,
                        });
                    let include_types: Vec<ReferenceType> =
                        include_types.into_iter().map(ReferenceType::from).collect();
                    let references = match alphabet {
                        IndexAlphabet::Dna => load_amrfinder_references(&db_dir, &include_types)?,
                        IndexAlphabet::Protein => {
                            load_amrfinder_protein_references(&db_dir, &include_types)?
                        }
                    };
                    let index = build_index(
                        &references,
                        &IndexBuildConfig {
                            alphabet,
                            k,
                            min_exact_gene_kmers,
                            min_hierarchy_unit_kmers,
                        },
                    )?;
                    save_index(&index, &out)?;
                    println!("{}", index.stats_string());
                    Ok(())
                }
                IndexCommand::Stats { index } => {
                    let index = load_index(&index)?;
                    print!("{}", index.stats_string());
                    Ok(())
                }
                IndexCommand::ReportMap { index, out } => {
                    let index = load_index(&index)?;
                    let text = report_map_tsv(&index);
                    if let Some(out) = out {
                        fs::write(&out, text)
                            .with_context(|| format!("write {}", out.display()))?;
                    } else {
                        print!("{text}");
                    }
                    Ok(())
                }
                IndexCommand::UnitStats { index, db_dir, out } => {
                    let index = load_index(&index)?;
                    let text = unit_stats_tsv(&index, db_dir.as_deref())?;
                    if let Some(out) = out {
                        fs::write(&out, text)
                            .with_context(|| format!("write {}", out.display()))?;
                    } else {
                        print!("{text}");
                    }
                    Ok(())
                }
            },
            CommandKind::Detect { command } => match command {
                DetectCommand::Direct {
                    index,
                    fasta,
                    sample_name,
                    params,
                } => run_detect(&index, &fasta, sample_name, QueryKind::Direct, params),
                DetectCommand::Cds {
                    index,
                    assembly,
                    cds_fasta,
                    protein,
                    protein_fasta,
                    out_dir,
                    orphos,
                    sample_name,
                    params,
                } => {
                    let sample = sample_name
                        .or_else(|| assembly.as_deref().map(sample_name_from_path))
                        .or_else(|| cds_fasta.as_deref().map(sample_name_from_path))
                        .or_else(|| protein_fasta.as_deref().map(sample_name_from_path))
                        .unwrap_or_else(|| "sample".to_string());
                    let input = if protein {
                        if let Some(protein_fasta) = protein_fasta {
                            protein_fasta
                        } else {
                            if cds_fasta.is_some() {
                                bail!(
                                    "detect cds --protein requires --assembly or --protein-fasta, not --cds-fasta"
                                );
                            }
                            let Some(assembly) = assembly else {
                                bail!(
                                    "detect cds --protein requires either --assembly or --protein-fasta"
                                );
                            };
                            let output = run_gene_caller(
                                &assembly,
                                &sample,
                                &GeneCallerConfig {
                                    out_dir,
                                    metagenomic: orphos.orphos_metagenomic,
                                    closed_ends: orphos.orphos_closed_ends,
                                    mask_n_runs: orphos.orphos_mask_n_runs,
                                    force_non_sd: orphos.orphos_force_non_sd,
                                    translation_table: orphos.orphos_translation_table,
                                },
                            )?;
                            output.protein_fasta
                        }
                    } else if let Some(cds) = cds_fasta {
                        cds
                    } else {
                        let Some(assembly) = assembly else {
                            bail!("detect cds requires either --cds-fasta or --assembly");
                        };
                        let output = run_gene_caller(
                            &assembly,
                            &sample,
                            &GeneCallerConfig {
                                out_dir,
                                metagenomic: orphos.orphos_metagenomic,
                                closed_ends: orphos.orphos_closed_ends,
                                mask_n_runs: orphos.orphos_mask_n_runs,
                                force_non_sd: orphos.orphos_force_non_sd,
                                translation_table: orphos.orphos_translation_table,
                            },
                        )?;
                        output.cds_fasta
                    };
                    let query_kind = if protein {
                        QueryKind::ProteinCds
                    } else {
                        QueryKind::Cds
                    };
                    run_detect(&index, &input, Some(sample), query_kind, params)
                }
            },
            CommandKind::Genes { command } => match command {
                GenesCommand::Call {
                    assembly,
                    out_dir,
                    orphos,
                    sample_name,
                } => {
                    let sample = sample_name.unwrap_or_else(|| sample_name_from_path(&assembly));
                    let output = run_gene_caller(
                        &assembly,
                        &sample,
                        &GeneCallerConfig {
                            out_dir,
                            metagenomic: orphos.orphos_metagenomic,
                            closed_ends: orphos.orphos_closed_ends,
                            mask_n_runs: orphos.orphos_mask_n_runs,
                            force_non_sd: orphos.orphos_force_non_sd,
                            translation_table: orphos.orphos_translation_table,
                        },
                    )?;
                    println!("{}", serde_json::to_string_pretty(&output)?);
                    Ok(())
                }
            },
            CommandKind::Eval { command } => match command {
                EvalCommand::NativeAmrfinder {
                    amrfinder_path,
                    db_dir,
                    assembly,
                    out,
                } => run_native_amrfinder(
                    amrfinder_path.as_deref().unwrap_or(Path::new("amrfinder")),
                    &db_dir,
                    &assembly,
                    &out,
                ),
                EvalCommand::DebugAmrfinderMisses {
                    index,
                    assembly,
                    amrfinder_tsv,
                    detector_json,
                    db_dir,
                    refinement_k,
                    missing_kmer_limit,
                    out,
                } => run_debug_amrfinder_misses(
                    &index,
                    &assembly,
                    &amrfinder_tsv,
                    &detector_json,
                    db_dir.as_deref(),
                    refinement_k,
                    missing_kmer_limit,
                    out.as_deref(),
                ),
                EvalCommand::TruthKmerEvidence {
                    index,
                    assembly,
                    amrfinder_tsv,
                    detector_json,
                    db_dir: _,
                    include_types,
                    min_gene_fraction,
                    min_family_fraction,
                    out,
                } => run_truth_kmer_evidence(
                    &index,
                    &assembly,
                    &amrfinder_tsv,
                    &detector_json,
                    &include_types,
                    min_gene_fraction,
                    min_family_fraction,
                    out.as_deref(),
                ),
            },
        }
    }

    fn report_map_tsv(index: &sparrowhawk_amr::AmrIndex) -> String {
        let mut out = String::from(
            "protein_accession\telement_symbol\ttype\tsubtype\thierarchy_node\treport_unit_key\treport_unit_type\treport_unit_id\treport_unit_label\n",
        );
        for gene in &index.genes {
            let unit = &index.units[gene.report_unit_id as usize];
            let unit_kind = unit.kind().as_str();
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}:{}\t{}\t{}\t{}\n",
                index.string(gene.protein_accession),
                index.string(gene.element_symbol),
                index.string(gene.type_name),
                index.string(gene.subtype),
                index.string(gene.hierarchy_node),
                unit_kind,
                index.string(unit.id),
                unit_kind,
                index.string(unit.id),
                index.string(unit.label),
            ));
        }
        out
    }

    #[derive(Debug, Default, Clone)]
    struct FamThresholds {
        complete_ident: String,
        complete_wp_coverage: String,
        complete_br_coverage: String,
        partial_ident: String,
        partial_wp_coverage: String,
        partial_br_coverage: String,
        reportable: String,
        family_name: String,
    }

    fn unit_stats_tsv(
        index: &sparrowhawk_amr::AmrIndex,
        db_dir: Option<&Path>,
    ) -> anyhow::Result<String> {
        let fam = if let Some(db_dir) = db_dir {
            load_fam_thresholds(&db_dir.join("fam.tsv"))?
        } else {
            std::collections::HashMap::new()
        };
        let mut out = String::from(
            "unit_id\tunit_key\tunit_type\tunit_label\telement_symbol\tgene_symbol\tallele_symbol\tgene_group\thierarchy_node\ttype\tsubtype\tclass\tsubclass\tdiagnostic_kmers\tmember_genes\tis_weak_hierarchy_unit\tblastrule_complete_ident\tblastrule_complete_wp_coverage\tblastrule_complete_br_coverage\tblastrule_partial_ident\tblastrule_partial_wp_coverage\tblastrule_partial_br_coverage\treportable\tfamily_name\n",
        );
        for (unit_idx, unit) in index.units.iter().enumerate() {
            let node = index.string(unit.hierarchy_node);
            let thresholds = fam.get(node).cloned().unwrap_or_default();
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                unit_idx,
                index.string(unit.id),
                unit.kind().as_str(),
                index.string(unit.label),
                index.optional_string(unit.element_symbol).unwrap_or_default(),
                index.optional_string(unit.gene_symbol).unwrap_or_default(),
                index.optional_string(unit.allele_symbol).unwrap_or_default(),
                index.string(unit.gene_group),
                node,
                index.string(unit.type_name),
                index.string(unit.subtype),
                index.string(unit.class_name),
                index.string(unit.subclass),
                unit.diagnostic_kmers,
                unit.member_count,
                unit.weak,
                thresholds.complete_ident,
                thresholds.complete_wp_coverage,
                thresholds.complete_br_coverage,
                thresholds.partial_ident,
                thresholds.partial_wp_coverage,
                thresholds.partial_br_coverage,
                thresholds.reportable,
                thresholds.family_name,
            ));
        }
        Ok(out)
    }

    fn load_fam_thresholds(
        path: &Path,
    ) -> anyhow::Result<std::collections::HashMap<String, FamThresholds>> {
        let text = fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
        let mut lines = text.lines();
        let Some(header) = lines.next() else {
            return Ok(std::collections::HashMap::new());
        };
        let columns: std::collections::HashMap<&str, usize> = header
            .trim_start_matches('#')
            .split('\t')
            .enumerate()
            .map(|(idx, name)| (name, idx))
            .collect();
        let mut by_node = std::collections::HashMap::new();
        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let node_id = fam_field(&fields, &columns, "node_id").to_string();
            if node_id.is_empty() {
                continue;
            }
            by_node.insert(
                node_id,
                FamThresholds {
                    complete_ident: fam_field(&fields, &columns, "blastrule_complete_ident")
                        .to_string(),
                    complete_wp_coverage: fam_field(
                        &fields,
                        &columns,
                        "blastrule_complete_wp_coverage",
                    )
                    .to_string(),
                    complete_br_coverage: fam_field(
                        &fields,
                        &columns,
                        "blastrule_complete_br_coverage",
                    )
                    .to_string(),
                    partial_ident: fam_field(&fields, &columns, "blastrule_partial_ident")
                        .to_string(),
                    partial_wp_coverage: fam_field(
                        &fields,
                        &columns,
                        "blastrule_partial_wp_coverage",
                    )
                    .to_string(),
                    partial_br_coverage: fam_field(
                        &fields,
                        &columns,
                        "blastrule_partial_br_coverage",
                    )
                    .to_string(),
                    reportable: fam_field(&fields, &columns, "reportable").to_string(),
                    family_name: fam_field(&fields, &columns, "family_name").to_string(),
                },
            );
        }
        Ok(by_node)
    }

    fn fam_field<'a>(
        fields: &'a [&str],
        columns: &std::collections::HashMap<&str, usize>,
        name: &str,
    ) -> &'a str {
        columns
            .get(name)
            .and_then(|idx| fields.get(*idx))
            .copied()
            .unwrap_or("")
    }

    fn fetch_db(out_dir: &Path, base_url: &str) -> anyhow::Result<()> {
        fs::create_dir_all(out_dir).with_context(|| format!("create {}", out_dir.display()))?;
        for file in DB_FILES {
            let url = format!("{}/{}", base_url.trim_end_matches('/'), file);
            let out = out_dir.join(file);
            let status = Command::new("curl")
                .arg("-L")
                .arg("-o")
                .arg(&out)
                .arg(&url)
                .status()
                .with_context(|| format!("run curl for {url}"))?;
            if !status.success() {
                bail!("curl failed for {url} with status {status}");
            }
        }
        Ok(())
    }

    fn run_detect(
        index_path: &Path,
        fasta_path: &Path,
        sample_name: Option<String>,
        query_kind: QueryKind,
        params: DetectArgs,
    ) -> anyhow::Result<()> {
        let index = load_index(index_path)?;
        let bytes =
            fs::read(fasta_path).with_context(|| format!("read {}", fasta_path.display()))?;
        let sample = sample_name.unwrap_or_else(|| sample_name_from_path(fasta_path));
        let params: DetectParams = params.into();
        let result = if query_kind == QueryKind::ProteinCds {
            detect_protein_fasta(&index, &bytes, &sample, &params)?
        } else {
            detect_fasta(&index, &bytes, &sample, query_kind, &params)?
        };
        println!("{}", serde_json::to_string_pretty(&result)?);
        Ok(())
    }

    fn run_native_amrfinder(
        amrfinder_path: &Path,
        db_dir: &Path,
        assembly: &Path,
        out: &Path,
    ) -> anyhow::Result<()> {
        let output = Command::new(amrfinder_path)
            .arg("-n")
            .arg(assembly)
            .arg("-d")
            .arg(db_dir)
            .output()
            .with_context(|| format!("run {}", amrfinder_path.display()))?;
        if !output.status.success() {
            bail!(
                "amrfinder failed with status {}: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr)
            );
        }
        fs::write(out, output.stdout).with_context(|| format!("write {}", out.display()))?;
        Ok(())
    }

    fn run_debug_amrfinder_misses(
        index_path: &Path,
        assembly: &Path,
        amrfinder_tsv: &Path,
        detector_json: &Path,
        db_dir: Option<&Path>,
        refinement_k: usize,
        missing_kmer_limit: usize,
        out: Option<&Path>,
    ) -> anyhow::Result<()> {
        let index = load_index(index_path)?;
        let report = debug_amrfinder_misses(DebugMissesConfig {
            index: &index,
            assembly_path: assembly,
            amrfinder_tsv,
            detector_json,
            db_dir,
            refinement_k,
            missing_kmer_limit,
        })?;
        let json = serde_json::to_string_pretty(&report)?;
        if let Some(out) = out {
            fs::write(out, json).with_context(|| format!("write {}", out.display()))?;
        } else {
            println!("{json}");
        }
        Ok(())
    }

    fn run_truth_kmer_evidence(
        index_path: &Path,
        assembly: &Path,
        amrfinder_tsv: &Path,
        detector_json: &Path,
        include_types: &[ReferenceTypeArg],
        min_gene_fraction: f64,
        min_family_fraction: f64,
        out: Option<&Path>,
    ) -> anyhow::Result<()> {
        let index = load_index(index_path)?;
        let include_types: Vec<ReferenceType> = include_types
            .iter()
            .copied()
            .map(ReferenceType::from)
            .collect();
        let report = sparrowhawk_amr::truth_kmer_evidence(TruthKmerEvidenceConfig {
            index: &index,
            assembly_path: assembly,
            amrfinder_tsv,
            detector_json,
            include_types: &include_types,
            min_gene_fraction,
            min_family_fraction,
        })?;
        let json = serde_json::to_string_pretty(&report)?;
        if let Some(out) = out {
            fs::write(out, json).with_context(|| format!("write {}", out.display()))?;
        } else {
            println!("{json}");
        }
        Ok(())
    }

    fn sample_name_from_path(path: &Path) -> String {
        path.file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("sample")
            .to_string()
    }
}

#[cfg(not(target_family = "wasm"))]
fn main() -> anyhow::Result<()> {
    native::run()
}
