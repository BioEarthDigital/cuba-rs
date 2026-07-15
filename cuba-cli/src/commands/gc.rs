use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct GcArgs {
    /// Input FASTA file
    pub fasta: String,
}

pub fn run(
    args: GcArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cuba_core::sequence::count_codons(&cds, &ct);
    let gc = cuba_core::metrics::gc::get_all_gc(&cf, &ct);

    let headers: Vec<String> = vec!["gene_id", "gc", "gc3s", "gc4d"]
        .iter().map(|s| s.to_string()).collect();
    let rows: Vec<Vec<String>> = gc
        .iter()
        .map(|r| {
            vec![
                r.gene_id.clone(),
                format!("{:.6}", r.gc),
                format!("{:.6}", r.gc3s),
                format!("{:.6}", r.gc4d),
            ]
        })
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
