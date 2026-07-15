use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct CaiArgs {
    /// Input FASTA file (target genes)
    pub fasta: String,

    /// Reference FASTA file (highly expressed genes for RSCU reference)
    #[arg(short = 'r', long)]
    pub reference: String,

    /// Analysis level: subfam (default) or amino_acid
    #[arg(short = 'l', long, default_value = "subfam")]
    pub level: String,
}

pub fn run(
    args: CaiArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;

    // Load reference genes and compute RSCU
    let ref_cds = super::load_cds(&args.reference, &ct)?;
    let ref_cf = cubar_core::sequence::count_codons(&ref_cds, &ct);
    let rscu = cubar_core::metrics::rscu::est_rscu(
        &ref_cf, None, 1.0, &ct, &args.level, false,
    );

    // Load target genes and compute CAI
    let target_cds = super::load_cds(&args.fasta, &ct)?;
    let target_cf = cubar_core::sequence::count_codons(&target_cds, &ct);
    let cai = cubar_core::metrics::cai::get_cai(&target_cf, &rscu, &args.level);

    let headers = vec!["gene_id".to_string(), "cai".to_string()];
    let rows: Vec<Vec<String>> = target_cf
        .gene_ids
        .iter()
        .enumerate()
        .map(|(i, id)| vec![id.clone(), format!("{:.6}", cai[i])])
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
