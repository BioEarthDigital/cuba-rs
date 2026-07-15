use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct RscuArgs {
    /// Input FASTA file
    pub fasta: String,

    /// Pseudo count (default: 1.0)
    #[arg(short = 'p', long, default_value = "1.0")]
    pub pseudo_cnt: f64,

    /// Analysis level: subfam (default) or amino_acid
    #[arg(short = 'l', long, default_value = "subfam")]
    pub level: String,

    /// Include stop codon RSCU values
    #[arg(long)]
    pub incl_stop: bool,
}

pub fn run(
    args: RscuArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cuba_core::sequence::count_codons(&cds, &ct);
    let rscu = cuba_core::metrics::rscu::est_rscu(
        &cf, None, args.pseudo_cnt, &ct, &args.level, args.incl_stop,
    );

    let headers: Vec<String> = vec![
        "amino_acid", "aa_code", "codon", "subfam", "count", "prop", "w_cai", "rscu",
    ].iter().map(|s| s.to_string()).collect();
    let rows: Vec<Vec<String>> = rscu
        .rows
        .iter()
        .map(|r| {
            vec![
                r.amino_acid.clone(),
                r.aa_code.to_string(),
                r.codon.clone(),
                r.subfam.clone(),
                format!("{:.2}", r.count),
                format!("{:.6}", r.prop),
                format!("{:.6}", r.w_cai),
                format!("{:.6}", r.rscu),
            ]
        })
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
