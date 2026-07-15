use std::collections::HashMap;
use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct TaiArgs {
    /// Input FASTA file
    pub fasta: String,

    /// tRNA gene copy number file (TSV: anticodon\tcopy_number)
    #[arg(short = 't', long)]
    pub trna_file: String,
}

pub fn run(
    args: TaiArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cuba_core::sequence::count_codons(&cds, &ct);

    // Read tRNA weights from file
    let content = std::fs::read_to_string(&args.trna_file)?;
    let mut trna_weights: HashMap<String, f64> = HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let anticodon = parts[0].trim().to_uppercase();
            let count: f64 = parts[1].trim().parse().unwrap_or(0.0);
            trna_weights.insert(anticodon, count);
        }
    }

    let tai = cuba_core::metrics::tai::get_tai(&cf, &trna_weights, &ct, None);

    let headers = vec!["gene_id".to_string(), "tai".to_string()];
    let rows: Vec<Vec<String>> = cf
        .gene_ids
        .iter()
        .enumerate()
        .map(|(i, id)| vec![id.clone(), format!("{:.6}", tai[i])])
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
