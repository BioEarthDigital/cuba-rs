use std::collections::HashMap;
use rayon::prelude::*;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

/// Calculate Fraction of Optimal Codons (Fop) for each gene.
///
/// Fop = number_of_optimal_codons / total_synonymous_codons
///
/// Optimal codons are defined by a list (typically derived from highly expressed genes).
///
/// # Arguments
/// * `cf` - Codon frequency matrix
/// * `optimal_codons` - Set of codons designated as "optimal"
/// * `codon_table` - Genetic code table
pub fn get_fop(
    cf: &CodonFreqMatrix,
    optimal_codons: &HashMap<String, bool>,
    codon_table: &CodonTable,
) -> Vec<f64> {
    cf.matrix
        .par_iter()
        .map(|row| {
            compute_fop_gene(row, &cf.codons, optimal_codons, codon_table)
        })
        .collect()
}

fn compute_fop_gene(
    row: &[f64],
    codon_names: &[String],
    optimal_codons: &HashMap<String, bool>,
    codon_table: &CodonTable,
) -> f64 {
    let mut optimal_count = 0.0f64;
    let mut total_syn_count = 0.0f64;

    for (j, &count) in row.iter().enumerate() {
        if count > 0.0 {
            let codon = &codon_names[j];

            // Skip stop codons
            if codon_table.is_stop(codon) {
                continue;
            }

            total_syn_count += count;

            if optimal_codons.contains_key(codon) {
                optimal_count += count;
            }
        }
    }

    if total_syn_count > 0.0 {
        optimal_count / total_syn_count
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};

    #[test]
    fn test_fop() {
        let ct = CodonTable::standard();

        // Gene: TTT (Phe), TTC (Phe, optimal), TTA (Leu), TTG (Leu, optimal)
        let cds = CdsSeq {
            id: "test".into(),
            seq: b"TTTTTCTTATTG".to_vec(),
            codons: vec!["TTT".into(), "TTC".into(), "TTA".into(), "TTG".into()],
        };

        let mut optimal = HashMap::new();
        optimal.insert("TTC".to_string(), true);
        optimal.insert("TTG".to_string(), true);

        let cf = count_codons(&[cds], &ct);
        let fop = get_fop(&cf, &optimal, &ct);

        // 2 optimal out of 4 synonymous = 0.5
        assert!((fop[0] - 0.5).abs() < 0.01, "Fop = {:.4}", fop[0]);
    }
}
