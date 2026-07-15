use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

pub use crate::sequence::{count_codons, count_codons_parallel};

/// Sum codons across all genes to get genome-wide counts
pub fn genome_wide_counts(cf: &CodonFreqMatrix, _codon_table: &CodonTable) -> Vec<(String, f64)> {
    let totals = cf.codon_totals();
    cf.codons
        .iter()
        .enumerate()
        .map(|(i, codon)| (codon.clone(), totals[i]))
        .collect()
}

/// Get total codon count per amino acid across all genes
pub fn amino_acid_counts(
    cf: &CodonFreqMatrix,
    codon_table: &CodonTable,
) -> Vec<(char, String, f64)> {
    let totals = cf.codon_totals();
    let mut aa_counts: std::collections::HashMap<char, (String, f64)> = std::collections::HashMap::new();

    for (i, codon) in cf.codons.iter().enumerate() {
        if let Some(info) = codon_table.codon_map.get(codon) {
            let entry = aa_counts.entry(info.aa_code).or_insert_with(|| {
                (info.amino_acid.clone(), 0.0)
            });
            entry.1 += totals[i];
        }
    }

    let mut result: Vec<_> = aa_counts
        .into_iter()
        .map(|(aa, (name, count))| (aa, name, count))
        .collect();
    result.sort_by_key(|(aa, _, _)| *aa);
    result
}
