use std::io::Write;

#[derive(clap::Args, Clone)]
pub struct OptimalArgs {
    /// Input FASTA file
    pub fasta: String,

    /// Optional expression data file (TSV: gene_id\texpression_level)
    #[arg(short = 'e', long)]
    pub expression: Option<String>,
}

pub fn run(
    args: OptimalArgs,
    gcid: &str,
    format: &str,
    writer: Box<dyn Write>,
) -> anyhow::Result<()> {
    let ct = super::load_codon_table(gcid)?;
    let cds = super::load_cds(&args.fasta, &ct)?;
    let cf = cuba_core::sequence::count_codons(&cds, &ct);

    let optimal = if let Some(expr_file) = &args.expression {
        // Read expression data
        let content = std::fs::read_to_string(expr_file)?;
        let mut expr_map: std::collections::HashMap<String, f64> = std::collections::HashMap::new();
        for line in content.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                let gene_id = parts[0].trim().to_string();
                let expr: f64 = parts[1].trim().parse().unwrap_or(0.0);
                expr_map.insert(gene_id, expr);
            }
        }

        // Build expression vector matching gene order in cf
        let expression: Vec<f64> = cf
            .gene_ids
            .iter()
            .map(|id| expr_map.get(id).copied().unwrap_or(0.0))
            .collect();

        cuba_core::optimize::est_optimal_codons_with_expression(&cf, &expression, &ct)
    } else {
        cuba_core::optimize::est_optimal_codons(&cf, &ct)
    };

    let headers: Vec<String> = vec!["amino_acid", "aa_code", "optimal_codon"]
        .iter().map(|s| s.to_string()).collect();
    let mut rows = Vec::new();
    for (aa, codon) in &optimal.optimal {
        rows.push(vec![
            cuba_core::genetic_code::aa1_to_aa3(*aa).to_string(),
            aa.to_string(),
            codon.clone(),
        ]);
    }
    rows.sort_by(|a, b| a[1].cmp(&b[1]));

    super::write_results(writer, format, &headers, &rows)
}
