use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct CountArgs {
    /// Input FASTA file containing coding sequences
    pub fasta: String,

    /// Group by amino acid instead of individual codons
    #[arg(long)]
    pub by_aa: bool,
}

pub fn run(
    args: CountArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cubar_core::sequence::count_codons(&cds, &ct);

    let mut headers: Vec<String> = vec!["gene_id".to_string()];
    headers.extend(cf.codons.iter().cloned());

    let rows: Vec<Vec<String>> = cf
        .gene_ids
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let mut row = vec![id.clone()];
            for count in &cf.matrix[i] {
                row.push(format!("{count}"));
            }
            row
        })
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
