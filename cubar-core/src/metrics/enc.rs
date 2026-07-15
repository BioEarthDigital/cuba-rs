use std::collections::HashMap;
use rayon::prelude::*;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

/// Calculate Effective Number of Codons (ENC) for each gene.
///
/// Uses the improved algorithm from Sun, Yang & Xia (2013) Mol Biol Evol 30:191-196,
/// as implemented in the cubar R package.
///
/// ENC ranges from 20 (extreme bias) to 61 (no bias).
///
/// Algorithm (matching cubar R):
/// 1. Group codons by subfamily (or amino_acid)
/// 2. For each gene and each group:
///    - n = sum of codon counts in this group
///    - p = (count + 1) / (n + k)  where k = number of codons in group
///    - f = sum(p^2)
/// 3. For each degeneracy class d:
///    - N_d = (n_groups_d) * sum(n) / sum(n * f)
///    - If sum(n) == 0, N_d = n_groups_d * d
/// 4. ENC = N_1 + N_2 + N_3 + N_4 + N_6 (includes non-degenerate = 1.0 each)
pub fn get_enc(cf: &CodonFreqMatrix, codon_table: &CodonTable, level: &str) -> Vec<f64> {
    let groups = codon_table.subfamily_groups();

    // Build list of codon groups (subfamilies), each as Vec<codon>
    let codon_groups: Vec<Vec<String>> = groups.into_values().collect();

    // Map codon -> group index
    let codon_to_group: HashMap<&str, usize> = {
        let mut m = HashMap::new();
        for (gi, codons) in codon_groups.iter().enumerate() {
            for c in codons {
                if !codon_table.is_stop(c) {
                    m.insert(c.as_str(), gi);
                }
            }
        }
        m
    };

    // Group degeneracies: number of codons in each group (excluding stops)
    let group_degs: Vec<usize> = codon_groups
        .iter()
        .map(|codons| codons.iter().filter(|c| !codon_table.is_stop(c)).count())
        .filter(|&d| d > 0)
        .collect();

    // Codon index in CF matrix
    let codon_to_idx: HashMap<&str, usize> = cf
        .codons
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();

    cf.matrix
        .par_iter()
        .map(|row| {
            compute_enc_gene(row, &codon_to_idx, &codon_to_group, &codon_groups, &group_degs)
        })
        .collect()
}

fn compute_enc_gene(
    row: &[f64],
    codon_to_idx: &HashMap<&str, usize>,
    codon_to_group: &HashMap<&str, usize>,
    codon_groups: &[Vec<String>],
    group_degs: &[usize],
) -> f64 {
    let n_groups = codon_groups.len();

    // For each group: compute n (total count) and f (sum of squared corrected frequencies)
    let mut n_per_group = vec![0.0f64; n_groups];
    let mut f_per_group = vec![0.0f64; n_groups];

    for (gi, codons) in codon_groups.iter().enumerate() {
        let k = codons.len() as f64;

        // Collect counts for codons in this group
        let counts: Vec<f64> = codons
            .iter()
            .map(|c| {
                codon_to_idx
                    .get(c.as_str())
                    .map(|&idx| row[idx])
                    .unwrap_or(0.0)
            })
            .collect();

        let n: f64 = counts.iter().sum();
        n_per_group[gi] = n;

        // p = (count + 1) / (n + k)
        // f = sum(p^2)
        let f: f64 = counts
            .iter()
            .map(|&cnt| {
                let p = (cnt + 1.0) / (n + k);
                p * p
            })
            .sum();
        f_per_group[gi] = f;
    }

    // Group by degeneracy (number of codons in group)
    let mut deg_to_groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for (gi, &deg) in group_degs.iter().enumerate() {
        deg_to_groups.entry(deg).or_default().push(gi);
    }

    // Collect unique degeneracies (sorted)
    let mut unique_degs: Vec<usize> = deg_to_groups.keys().copied().collect();
    unique_degs.sort();

    // ENC = sum over degeneracy classes
    let mut enc = 0.0;

    for &deg in &unique_degs {
        let group_indices = &deg_to_groups[&deg];
        let n_groups_deg = group_indices.len();

        let sum_n: f64 = group_indices.iter().map(|&gi| n_per_group[gi]).sum();
        let sum_nf: f64 = group_indices.iter().map(|&gi| n_per_group[gi] * f_per_group[gi]).sum();

        if sum_n == 0.0 {
            // If no codons for this degeneracy class, use theoretical max
            enc += (n_groups_deg * deg) as f64;
        } else {
            // N_deg = n_groups_deg * sum(n) / sum(n * f)
            enc += (n_groups_deg as f64) * sum_n / sum_nf;
        }
    }

    enc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};

    #[test]
    fn test_enc_uniform() {
        let ct = CodonTable::standard();

        // Create a gene that uses all synonymous codons equally — use many repeats
        let groups = ct.subfamily_groups();
        let mut seq = String::new();
        for (_sf, codons) in &groups {
            for _ in 0..10 {
                for codon in codons {
                    seq.push_str(codon);
                }
            }
        }

        let cds = CdsSeq {
            id: "uniform".into(),
            seq: seq.as_bytes().to_vec(),
            codons: seq.as_bytes().chunks(3).map(|c| String::from_utf8_lossy(c).to_string()).collect(),
        };

        let cf = count_codons(&[cds], &ct);
        let enc = get_enc(&cf, &ct, "subfam");

        // With near-uniform usage, ENC should be high
        assert!(enc[0] > 48.0, "Expected ENC > 48 for uniform usage, got {}", enc[0]);
    }

    #[test]
    fn test_enc_extreme_bias() {
        let ct = CodonTable::standard();

        // Gene using only one codon per amino acid subfamily, repeated many times
        let groups = ct.subfamily_groups();
        let mut seq = String::new();
        for (_sf, codons) in &groups {
            if !codons.is_empty() {
                for _ in 0..10 {
                    seq.push_str(&codons[0]);
                }
            }
        }

        let cds = CdsSeq {
            id: "biased".into(),
            seq: seq.as_bytes().to_vec(),
            codons: seq.as_bytes().chunks(3).map(|c| String::from_utf8_lossy(c).to_string()).collect(),
        };

        let cf = count_codons(&[cds], &ct);
        let enc = get_enc(&cf, &ct, "subfam");

        // With extreme bias, ENC should be lower
        assert!(enc[0] < 40.0, "Expected ENC < 40 for extreme bias, got {}", enc[0]);
    }

    #[test]
    fn test_enc_single_gene() {
        let ct = CodonTable::standard();

        // A simple gene: ATG GCT GGT TAA
        let cds = CdsSeq {
            id: "simple".into(),
            seq: b"ATGGCTGGTTAA".to_vec(),
            codons: vec!["ATG".into(), "GCT".into(), "GGT".into(), "TAA".into()],
        };

        let cf = count_codons(&[cds], &ct);
        let enc = get_enc(&cf, &ct, "subfam");

        // ENC should be a reasonable value
        assert!(enc[0] >= 0.0, "ENC should be non-negative");
        assert!(enc[0] <= 61.0, "ENC should be <= 61");
    }
}
