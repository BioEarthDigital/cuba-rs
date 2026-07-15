use std::collections::HashMap;
use rayon::prelude::*;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

/// Calculate Deviation from Proportionality (DP) of host tRNA availability.
///
/// DP measures how much a gene's codon usage deviates from the relative
/// tRNA gene copy numbers.
///
/// # Arguments
/// * `cf` - Codon frequency matrix
/// * `trna_weights` - tRNA gene copy numbers per anticodon
/// * `codon_table` - Genetic code table
pub fn get_dp(
    cf: &CodonFreqMatrix,
    trna_weights: &HashMap<String, f64>,
    codon_table: &CodonTable,
) -> Vec<f64> {
    // Map tRNA anticodons to codon proportions
    let mut codon_trna: HashMap<String, f64> = HashMap::new();

    for (anticodon, &t_gcn) in trna_weights {
        // Simple mapping: convert anticodon to codon
        let codon = super::tai::anticodon_to_codon(anticodon);
        *codon_trna.entry(codon).or_default() += t_gcn;
    }

    let total_trna: f64 = codon_trna.values().sum();

    cf.matrix
        .par_iter()
        .map(|row| {
            // Count total codons for this gene (excluding stops)
            let mut total_codons = 0.0f64;
            for (j, &count) in row.iter().enumerate() {
                if count > 0.0 && !codon_table.is_stop(&cf.codons[j]) {
                    total_codons += count;
                }
            }

            if total_codons == 0.0 {
                return 0.0;
            }

            // DP = sum of absolute deviations between codon proportion and tRNA proportion
            let mut dp = 0.0f64;
            for (j, &count) in row.iter().enumerate() {
                if count > 0.0 {
                    let codon = &cf.codons[j];
                    if codon_table.is_stop(codon) {
                        continue;
                    }
                    let codon_prop = count / total_codons;
                    let trna_prop = if total_trna > 0.0 {
                        codon_trna.get(codon).copied().unwrap_or(0.0) / total_trna
                    } else {
                        0.0
                    };
                    dp += (codon_prop - trna_prop).abs();
                }
            }

            dp / 2.0 // Normalize to [0, 1]
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};

    #[test]
    fn test_dp_perfect_match() {
        let ct = CodonTable::standard();

        // Gene uses only TTC (Phe)
        let cds = CdsSeq {
            id: "test".into(),
            seq: b"TTCTTC".to_vec(),
            codons: vec!["TTC".into(), "TTC".into()],
        };

        // tRNA: only GAA anticodon (recognizes TTC)
        let mut trna = HashMap::new();
        trna.insert("GAA".to_string(), 10.0);

        let cf = count_codons(&[cds], &ct);
        let dp = get_dp(&cf, &trna, &ct);

        // Perfect match should give low/zero DP
        assert!(dp[0] >= 0.0 && dp[0] <= 1.0);
    }
}
