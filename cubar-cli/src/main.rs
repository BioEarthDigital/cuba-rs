mod commands;

use clap::{Parser, Subcommand};
use commands::count::CountArgs;
use commands::enc::EncArgs;
use commands::cai::CaiArgs;
use commands::rscu::RscuArgs;
use commands::fop::FopArgs;
use commands::tai::TaiArgs;
use commands::gc::GcArgs;
use commands::optimal::OptimalArgs;
use commands::optimize::OptimizeArgs;
use commands::slide::SlideArgs;
use commands::show_code::ShowCodeArgs;

/// cubar — Codon Usage Bias Analysis in Rust
///
/// A fast, memory-efficient tool for analyzing codon usage bias in coding sequences.
#[derive(Parser)]
#[command(name = "cubar", version, about, long_about = None)]
struct Cli {
    /// NCBI genetic code ID (default: 1 = Standard)
    #[arg(short = 'c', long, global = true, default_value = "1")]
    gcid: String,

    /// Output format: csv, tsv, or json
    #[arg(short = 'f', long, global = true, default_value = "csv")]
    format: String,

    /// Output file (default: stdout)
    #[arg(short = 'o', long, global = true)]
    output: Option<String>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Count codon frequencies in coding sequences
    Count(CountArgs),

    /// Calculate Effective Number of Codons (ENC)
    Enc(EncArgs),

    /// Calculate Codon Adaptation Index (CAI)
    Cai(CaiArgs),

    /// Estimate Relative Synonymous Codon Usage (RSCU)
    Rscu(RscuArgs),

    /// Calculate Fraction of Optimal Codons (Fop)
    Fop(FopArgs),

    /// Calculate tRNA Adaptation Index (tAI)
    Tai(TaiArgs),

    /// Calculate GC content metrics (GC, GC3s, GC4d)
    Gc(GcArgs),

    /// Identify optimal codons
    Optimal(OptimalArgs),

    /// Optimize codon usage of a sequence
    Optimize(OptimizeArgs),

    /// Sliding window analysis
    Slide(SlideArgs),

    /// List available genetic codes
    ListCodes,

    /// Show genetic code table
    ShowCode(ShowCodeArgs),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    // Set up output
    let output: Box<dyn std::io::Write> = if let Some(path) = &cli.output {
        let file = std::fs::File::create(path)?;
        Box::new(std::io::BufWriter::new(file))
    } else {
        Box::new(std::io::stdout().lock())
    };

    match cli.command {
        Commands::Count(args) => commands::count::run(args, &cli.gcid, &cli.format, output),
        Commands::Enc(args) => commands::enc::run(args, &cli.gcid, &cli.format, output),
        Commands::Cai(args) => commands::cai::run(args, &cli.gcid, &cli.format, output),
        Commands::Rscu(args) => commands::rscu::run(args, &cli.gcid, &cli.format, output),
        Commands::Fop(args) => commands::fop::run(args, &cli.gcid, &cli.format, output),
        Commands::Tai(args) => commands::tai::run(args, &cli.gcid, &cli.format, output),
        Commands::Gc(args) => commands::gc::run(args, &cli.gcid, &cli.format, output),
        Commands::Optimal(args) => commands::optimal::run(args, &cli.gcid, &cli.format, output),
        Commands::Optimize(args) => commands::optimize::run(args, &cli.gcid, &cli.format, output),
        Commands::Slide(args) => commands::slide::run(args, &cli.gcid, &cli.format, output),
        Commands::ListCodes => commands::list_codes(),
        Commands::ShowCode(args) => commands::show_code::run(args),
    }
}
