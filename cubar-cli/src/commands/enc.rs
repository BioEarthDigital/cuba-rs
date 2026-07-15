use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct EncArgs {
    /// Input FASTA file
    pub fasta: String,

    /// Analysis level: subfam (default) or amino_acid
    #[arg(short = 'l', long, default_value = "subfam")]
    pub level: String,
}

pub fn run(
    args: EncArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cubar_core::sequence::count_codons(&cds, &ct);
    let enc = cubar_core::metrics::enc::get_enc(&cf, &ct, &args.level);

    let headers = vec!["gene_id".to_string(), "enc".to_string()];
    let rows: Vec<Vec<String>> = cf
        .gene_ids
        .iter()
        .enumerate()
        .map(|(i, id)| vec![id.clone(), format!("{:.4}", enc[i])])
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
