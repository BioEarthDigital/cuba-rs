use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct OptimizeArgs {
    /// Input FASTA file containing sequence(s) to optimize
    pub fasta: String,

    /// File listing optimal codons (one per line)
    #[arg(short = 'O', long)]
    pub optimal_codons: String,

    /// Optimization method: naive (default)
    #[arg(short = 'm', long, default_value = "naive")]
    pub method: String,
}

pub fn run(
    args: OptimizeArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let _cf = cuba_core::sequence::count_codons(&cds, &ct);

    // Read optimal codons from file
    let content = std::fs::read_to_string(&args.optimal_codons)?;
    let optimal_codons_list: Vec<String> = content
        .lines()
        .map(|l| l.trim().to_uppercase())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();

    // Build OptimalCodonTable from the list
    let mut optimal = std::collections::HashMap::new();
    let mut optimal_list: std::collections::HashMap<char, Vec<String>> = std::collections::HashMap::new();

    for codon in &optimal_codons_list {
        if let Some(info) = ct.codon_map.get(codon) {
            optimal.entry(info.aa_code).or_insert_with(|| codon.clone());
            optimal_list.entry(info.aa_code).or_default().push(codon.clone());
        }
    }

    let opt_table = cuba_core::optimize::OptimalCodonTable { optimal, optimal_list };

    let headers: Vec<String> = vec!["gene_id", "original_length", "optimized_seq"]
        .iter().map(|s| s.to_string()).collect();
    let mut rows = Vec::new();

    for seq in &cds {
        let original = String::from_utf8_lossy(&seq.seq).to_string();
        let optimized = cuba_core::optimize::codon_optimize(&original, &opt_table, &args.method);
        rows.push(vec![
            seq.id.clone(),
            original.len().to_string(),
            optimized,
        ]);
    }

    super::write_results(writer, format, &headers, &rows)
}
