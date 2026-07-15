use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct SlideArgs {
    /// Input FASTA file
    pub fasta: String,

    /// Window size in codons
    #[arg(short = 'w', long, default_value = "20")]
    pub window_size: usize,

    /// Step size in codons
    #[arg(short = 's', long, default_value = "1")]
    pub step_size: usize,

    /// Metrics to compute: enc, cai, gc, all
    #[arg(short = 'M', long, default_value = "all")]
    pub metric: String,

    /// Reference FASTA file for CAI calculation (required if metric includes cai)
    #[arg(short = 'r', long)]
    pub reference: Option<String>,
}

pub fn run(
    args: SlideArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;

    // Compute RSCU from reference if needed
    let rscu = if let Some(ref_path) = &args.reference {
        let ref_cds = super::load_cds(ref_path, &ct)?;
        let ref_cf = cuba_core::sequence::count_codons(&ref_cds, &ct);
        Some(cuba_core::metrics::rscu::est_rscu(
            &ref_cf, None, 1.0, &ct, "subfam", false,
        ))
    } else {
        None
    };

    let metrics: Vec<&str> = args.metric.split(',').collect();
    let results = cuba_core::slide::slide(
        &cds,
        args.window_size,
        args.step_size,
        &ct,
        rscu.as_ref(),
        &metrics,
    );

    let headers: Vec<String> = vec![
        "seq_id", "start", "end", "window_index", "enc", "cai", "gc", "gc3s", "gc4d",
    ].iter().map(|s| s.to_string()).collect();
    let rows: Vec<Vec<String>> = results
        .iter()
        .map(|w| {
            vec![
                w.seq_id.clone(),
                w.start.to_string(),
                w.end.to_string(),
                w.window_index.to_string(),
                w.enc.map(|v| format!("{:.4}", v)).unwrap_or_default(),
                w.cai.map(|v| format!("{:.6}", v)).unwrap_or_default(),
                w.gc.map(|v| format!("{:.6}", v)).unwrap_or_default(),
                w.gc3s.map(|v| format!("{:.6}", v)).unwrap_or_default(),
                w.gc4d.map(|v| format!("{:.6}", v)).unwrap_or_default(),
            ]
        })
        .collect();

    super::write_results(writer, format, &headers, &rows)
}
