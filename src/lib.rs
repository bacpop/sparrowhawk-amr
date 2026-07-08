pub mod amrfinder_db;
#[cfg(not(target_family = "wasm"))]
pub mod debug;
pub mod detect;
pub mod fasta;
#[cfg(not(target_family = "wasm"))]
pub mod gene_callers;
pub mod index;
pub mod kmer;
pub mod translate;

pub use amrfinder_db::{
    AmrReference, HierarchyNode, ReferenceType, load_amrfinder_protein_references,
    load_amrfinder_references,
};
#[cfg(not(target_family = "wasm"))]
pub use debug::{
    DebugMissesConfig, DebugMissesReport, TruthKmerEvidenceConfig, TruthKmerEvidenceReport,
    debug_amrfinder_misses, truth_kmer_evidence,
};
pub use detect::{
    DetectParams, DetectionResult, QueryKind, RefinementMode, detect_fasta, detect_protein_fasta,
};
#[cfg(not(target_family = "wasm"))]
pub use gene_callers::{GeneCaller, GeneCallerConfig, run_gene_caller};
pub use index::{
    AmrIndex, IndexAlphabet, IndexBuildConfig, ReportUnit, ReportUnitKind, build_index,
    load_index_from_bytes,
};
#[cfg(not(target_family = "wasm"))]
pub use index::{load_index, save_index};
pub use translate::{DEFAULT_BACTERIAL_TRANSLATION_TABLE, translate_cds};

#[cfg(target_family = "wasm")]
use wasm_bindgen::prelude::*;

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Main interface for webassembly
#[cfg(target_family = "wasm")]
#[wasm_bindgen]
pub struct AmrDetector {
    index: AmrIndex,
}

#[cfg(target_family = "wasm")]
#[wasm_bindgen]
impl AmrDetector {
    #[wasm_bindgen(constructor)]
    pub fn new(index_bytes: &[u8]) -> Result<AmrDetector, JsValue> {
        init_panic_hook();
        let index = load_index_from_bytes(index_bytes).map_err(js_error)?;
        if index.alphabet != IndexAlphabet::Dna {
            return Err(JsValue::from_str(
                "AMR direct mode requires a DNA AMR index",
            ));
        }
        Ok(Self { index })
    }

    pub fn info(&self) -> String {
        self.index.stats_string()
    }

    /// Detect directly on contigs, not on called genes.
    pub fn detect_direct(
        &self,
        sample_name: &str,
        fasta_bytes: &[u8],
        min_gene_fraction: f64,
        min_gene_group_fraction: f64,
    ) -> Result<String, JsValue> {
        validate_fraction("min_gene_fraction", min_gene_fraction)?; // Probably not needed...
        validate_fraction("min_gene_group_fraction", min_gene_group_fraction)?; // Probably not needed...
        let params = DetectParams {
            min_gene_fraction,
            min_gene_group_fraction,
            ..DetectParams::default()
        };
        let result = detect_fasta(
            &self.index,
            fasta_bytes,
            sample_name,
            QueryKind::Direct,
            &params,
        )
        .map_err(js_error)?;
        serde_json::to_string(&result).map_err(js_error)
    }

    pub fn detect_cds(
        &self,
        sample_name: &str,
        fasta_bytes: &[u8],
        min_gene_fraction: f64,
        min_gene_group_fraction: f64,
    ) -> Result<String, JsValue> {
        validate_fraction("min_gene_fraction", min_gene_fraction)?; // Probably not needed...
        validate_fraction("min_gene_group_fraction", min_gene_group_fraction)?; // Probably not needed...
        let params = DetectParams {
            min_gene_fraction,
            min_gene_group_fraction,
            ..DetectParams::default()
        };
        let result = detect_fasta(
            &self.index,
            fasta_bytes,
            sample_name,
            QueryKind::Cds,
            &params,
        )
        .map_err(js_error)?;
        serde_json::to_string(&result).map_err(js_error)
    }
}

#[cfg(target_family = "wasm")]
fn validate_fraction(name: &str, value: f64) -> Result<(), JsValue> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        return Ok(());
    }
    Err(JsValue::from_str(&format!(
        "{name} must be between 0 and 1"
    )))
}

#[cfg(target_family = "wasm")]
fn js_error(error: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&error.to_string())
}
