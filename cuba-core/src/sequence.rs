use std::collections::HashMap;
use std::path::Path;
use rayon::prelude::*;
use anyhow::{Context, Result, bail};

use crate::genetic_code::CodonTable;

/// A validated coding sequence
#[derive(Debug, Clone)]
pub struct CdsSeq {
    /// Sequence identifier (from FASTA header)
    pub id: String,
    /// Uppercase DNA sequence (without newlines/whitespace)
    pub seq: Vec<u8>,
    /// Codons (triplets)
    pub codons: Vec<String>,
}

impl CdsSeq {
    /// Number of codons (excluding stop)
    pub fn n_codons(&self) -> usize {
        self.codons.len()
    }

    /// GC content of the full CDS
    pub fn gc_content(&self) -> f64 {
        if self.seq.is_empty() {
            return 0.0;
        }
        let gc = self.seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc as f64 / self.seq.len() as f64
    }
}

/// Read FASTA file and return a vector of sequences (raw, not yet split into codons).
pub fn read_fasta(path: &Path) -> Result<Vec<CdsSeq>> {
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("Cannot open FASTA file: {}", path.display()))?;

    let mut seqs = Vec::new();

    while let Some(record) = reader.next() {
        let record = record.with_context(|| "Failed to read FASTA record")?;
        let id = String::from_utf8_lossy(record.id()).to_string();
        // Take the first word of the header as the ID
        let id = id.split_whitespace().next().unwrap_or(&id).to_string();
        let seq: Vec<u8> = record.seq().to_ascii_uppercase();

        seqs.push(CdsSeq {
            id,
            seq,
            codons: Vec::new(),
        });
    }

    Ok(seqs)
}

/// Validate coding sequences:
/// - Length must be a multiple of 3
/// - Must start with a start codon (warning if not)
/// - Must end with a stop codon (warning if not)
/// - Must not contain internal stop codons
///
/// Returns the validated sequences with codons populated.
pub fn check_cds(seqs: Vec<CdsSeq>, codon_table: &CodonTable) -> Result<Vec<CdsSeq>> {
    let mut valid = Vec::new();
    let mut skipped = 0u64;

    for mut s in seqs {
        if s.seq.len() % 3 != 0 {
            eprintln!("Warning: sequence {} length ({}) is not a multiple of 3, skipping", s.id, s.seq.len());
            skipped += 1;
            continue;
        }

        let codons: Vec<String> = seq_to_codons(&s.seq);
        if codons.is_empty() {
            skipped += 1;
            continue;
        }

        // Check for internal stop codons (all except the last codon)
        let mut has_internal_stop = false;
        for (i, codon) in codons.iter().enumerate() {
            let is_last = i == codons.len() - 1;
            if codon_table.is_stop(codon) && !is_last {
                eprintln!(
                    "Warning: sequence {} has internal stop codon {} at position {}",
                    s.id,
                    codon,
                    i * 3 + 1
                );
                has_internal_stop = true;
                break;
            }
        }

        if has_internal_stop {
            skipped += 1;
            continue;
        }

        // Warn about missing start codon
        if !codon_table.is_start(&codons[0]) {
            eprintln!(
                "Note: sequence {} does not start with a standard start codon (starts with {})",
                s.id, codons[0]
            );
        }

        // Warn about missing terminal stop codon
        let last = codons.last().unwrap();
        if !codon_table.is_stop(last) {
            eprintln!(
                "Note: sequence {} does not end with a stop codon (ends with {})",
                s.id, last
            );
        }

        s.codons = codons;
        valid.push(s);
    }

    if valid.is_empty() {
        bail!("No valid coding sequences found ({} skipped)", skipped);
    }

    eprintln!("Validated {} coding sequences ({} skipped)", valid.len(), skipped);
    Ok(valid)
}

/// Split a DNA sequence into codon triplets.
pub fn seq_to_codons(seq: &[u8]) -> Vec<String> {
    seq.chunks(3)
        .map(|chunk| {
            String::from_utf8_lossy(chunk).to_ascii_uppercase()
        })
        .collect()
}

/// A codon frequency matrix: rows = genes, columns = codons.
#[derive(Debug, Clone)]
pub struct CodonFreqMatrix {
    /// Gene IDs (row names)
    pub gene_ids: Vec<String>,
    /// Codon names (column names), in canonical order
    pub codons: Vec<String>,
    /// Matrix: matrix\[gene_idx\]\[codon_idx\] = count
    pub matrix: Vec<Vec<f64>>,
}

impl CodonFreqMatrix {
    /// Total count per gene
    pub fn gene_totals(&self) -> Vec<f64> {
        self.matrix
            .iter()
            .map(|row| row.iter().sum())
            .collect()
    }

    /// Total count per codon (across all genes)
    pub fn codon_totals(&self) -> Vec<f64> {
        let n_codons = self.codons.len();
        let mut totals = vec![0.0; n_codons];
        for row in &self.matrix {
            for (j, &count) in row.iter().enumerate() {
                totals[j] += count;
            }
        }
        totals
    }

    /// Number of genes (rows)
    pub fn n_genes(&self) -> usize {
        self.gene_ids.len()
    }

    /// Number of codons (columns)
    pub fn n_codons(&self) -> usize {
        self.codons.len()
    }
}

/// Count codon frequencies for a set of CDS sequences.
///
/// Returns a CodonFreqMatrix where each row corresponds to a gene
/// and each column corresponds to a codon.
pub fn count_codons(cds_seqs: &[CdsSeq], codon_table: &CodonTable) -> CodonFreqMatrix {
    // Do NOT exclude stop codons — the caller (e.g., est_rscu) decides
    let codons: Vec<String> = codon_table.all_codons.clone();
    let codon_to_idx: HashMap<&str, usize> = codons
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();

    let n_codons = codons.len();
    let mut matrix = Vec::with_capacity(cds_seqs.len());
    let mut gene_ids = Vec::with_capacity(cds_seqs.len());

    for seq in cds_seqs {
        let mut row = vec![0.0; n_codons];
        for codon in &seq.codons {
            if let Some(&idx) = codon_to_idx.get(codon.as_str()) {
                row[idx] += 1.0;
            }
        }
        matrix.push(row);
        gene_ids.push(seq.id.clone());
    }

    CodonFreqMatrix {
        gene_ids,
        codons,
        matrix,
    }
}

/// Count codon frequencies in parallel using rayon.
pub fn count_codons_parallel(cds_seqs: &[CdsSeq], codon_table: &CodonTable) -> CodonFreqMatrix {
    let codons: Vec<String> = codon_table.all_codons.clone();
    let codon_to_idx: HashMap<&str, usize> = codons
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();

    let n_codons = codons.len();

    let results: Vec<(String, Vec<f64>)> = cds_seqs
        .par_iter()
        .map(|seq| {
            let mut row = vec![0.0; n_codons];
            for codon in &seq.codons {
                if let Some(&idx) = codon_to_idx.get(codon.as_str()) {
                    row[idx] += 1.0;
                }
            }
            (seq.id.clone(), row)
        })
        .collect();

    let (gene_ids, matrix): (Vec<_>, Vec<_>) = results.into_iter().unzip();

    CodonFreqMatrix {
        gene_ids,
        codons,
        matrix,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;

    #[test]
    fn test_seq_to_codons() {
        let seq = b"ATGGGCTAA";
        let codons = seq_to_codons(seq);
        assert_eq!(codons, vec!["ATG", "GGC", "TAA"]);
    }

    #[test]
    fn test_count_codons() {
        let ct = CodonTable::standard();

        let cds = CdsSeq {
            id: "test".into(),
            seq: b"ATGTGGTAA".to_vec(),
            codons: vec!["ATG".into(), "TGG".into(), "TAA".into()],
        };

        let cf = count_codons(&[cds], &ct);

        let atg_idx = cf.codons.iter().position(|c| c == "ATG").unwrap();
        let tgg_idx = cf.codons.iter().position(|c| c == "TGG").unwrap();
        let taa_idx = cf.codons.iter().position(|c| c == "TAA").unwrap();

        assert_eq!(cf.matrix[0][atg_idx], 1.0);
        assert_eq!(cf.matrix[0][tgg_idx], 1.0);
        assert_eq!(cf.matrix[0][taa_idx], 1.0);

        // Other codons should be zero
        let total: f64 = cf.matrix[0].iter().sum();
        assert_eq!(total, 3.0);
    }

    #[test]
    fn test_gc_content() {
        let cds = CdsSeq {
            id: "test".into(),
            seq: b"GCGCATATG".to_vec(), // G,C,G,C,A,T,A,T,G = 9 bases, 5 GC
            codons: vec![],
        };
        assert!((cds.gc_content() - 5.0/9.0).abs() < 0.01);
    }
}
