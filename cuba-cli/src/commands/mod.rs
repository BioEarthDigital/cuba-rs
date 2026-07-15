pub mod count;
pub mod enc;
pub mod cai;
pub mod rscu;
pub mod fop;
pub mod tai;
pub mod gc;
pub mod optimal;
pub mod optimize;
pub mod slide;

use cuba_core::genetic_code::CodonTable;

/// Load a genetic code table by ID
pub fn load_codon_table(gcid: &str) -> anyhow::Result<CodonTable> {
    CodonTable::from_ncbi_id(gcid)
        .map_err(|e| anyhow::anyhow!("{e}"))
}

/// Read and validate coding sequences from a FASTA file
pub fn load_cds(path: &str, codon_table: &CodonTable) -> anyhow::Result<Vec<cuba_core::sequence::CdsSeq>> {
    let seqs = cuba_core::sequence::read_fasta(std::path::Path::new(path))?;
    cuba_core::sequence::check_cds(seqs, codon_table)
}

/// Write results in CSV, TSV, or JSON format
pub fn write_results(
    writer: Box<dyn std::io::Write>,
    format: &str,
    headers: &[String],
    rows: &[Vec<String>],
) -> anyhow::Result<()> {
    match format {
        "csv" => write_csv(writer, headers, rows),
        "tsv" => write_tsv(writer, headers, rows),
        "json" => write_json(writer, headers, rows),
        _ => anyhow::bail!("Unsupported format: {format}. Use csv, tsv, or json."),
    }
}

fn write_csv(
    mut writer: Box<dyn std::io::Write>,
    headers: &[String],
    rows: &[Vec<String>],
) -> anyhow::Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .has_headers(true)
        .from_writer(Vec::new());
    wtr.write_record(headers)?;
    for row in rows {
        wtr.write_record(row)?;
    }
    let data = wtr.into_inner()?;
    writer.write_all(&data)?;
    Ok(())
}

fn write_tsv(
    mut writer: Box<dyn std::io::Write>,
    headers: &[String],
    rows: &[Vec<String>],
) -> anyhow::Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(Vec::new());
    wtr.write_record(headers)?;
    for row in rows {
        wtr.write_record(row)?;
    }
    let data = wtr.into_inner()?;
    writer.write_all(&data)?;
    Ok(())
}

fn write_json(
    mut writer: Box<dyn std::io::Write>,
    headers: &[String],
    rows: &[Vec<String>],
) -> anyhow::Result<()> {
    let objects: Vec<serde_json::Value> = rows
        .iter()
        .map(|row| {
            let mut obj = serde_json::Map::new();
            for (h, val) in headers.iter().zip(row.iter()) {
                obj.insert(h.clone(), serde_json::Value::String(val.clone()));
            }
            serde_json::Value::Object(obj)
        })
        .collect();
    serde_json::to_writer_pretty(&mut writer, &objects)?;
    Ok(())
}

/// List available genetic codes
pub fn list_codes() -> anyhow::Result<()> {
    println!("{:<6} {}", "ID", "Name");
    println!("{}", "-".repeat(70));
    for (id, name) in CodonTable::list_available() {
        println!("{:<6} {name}", id);
    }
    Ok(())
}

/// Show a specific genetic code table
pub mod show_code {
    use cuba_core::genetic_code::CodonTable;

    #[derive(clap::Args, Clone)]
    pub struct ShowCodeArgs {
        /// NCBI genetic code ID
        #[arg(default_value = "1")]
        pub gcid: String,
    }

    pub fn run(args: ShowCodeArgs) -> anyhow::Result<()> {
        let ct = CodonTable::from_ncbi_id(&args.gcid)
            .map_err(|e| anyhow::anyhow!("{e}"))?;

        println!("Genetic Code: {} (ID: {})", ct.name, ct.gcid);
        println!();
        println!("{:<10} {:<8} {:<6} {:<20}", "Codon", "AA", "AA3", "Subfamily");
        println!("{}", "-".repeat(50));

        for codon in &ct.all_codons {
            if let Some(info) = ct.codon_map.get(codon) {
                let start = if ct.is_start(codon) { " (start)" } else { "" };
                let stop = if ct.is_stop(codon) { " (stop)" } else { "" };
                println!(
                    "{:<10} {:<8} {:<6} {:<20}{}{}",
                    codon, info.aa_code, info.amino_acid, info.subfam, start, stop
                );
            }
        }

        Ok(())
    }
}
