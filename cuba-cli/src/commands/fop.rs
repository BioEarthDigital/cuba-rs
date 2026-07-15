use std::collections::HashMap;
use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct FopArgs {
    /// Input FASTA file
    pub fasta: String,

    /// File listing optimal codons (one per line)
    #[arg(short = 'O', long)]
    pub optimal_codons: String,
}

pub fn run(
    args: FopArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cuba_core::sequence::count_codons(&cds, &ct);

    // Read optimal codons from file
    let content = std::fs::read_to_string(&args.optimal_codons)?;
    let optimal: HashMap<String, bool> = content
        .lines()
        .map(|l| l.trim().to_uppercase())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| (l, true))
        .collect();

    let fop = cuba_core::metrics::fop::get_fop(&cf, &optimal, &ct);

    let headers = vec!["gene_id".to_string(), "fop".to_string()];
    let rows: Vec<Vec<String>> = cf
        .gene_ids
        .iter()
        .enumerate()
        .map(|(i, id)| vec![id.clone(), format!("{:.6}", fop[i])])
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
