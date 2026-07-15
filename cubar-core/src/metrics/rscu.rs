use std::collections::HashMap;
use serde::{Deserialize, Serialize};

use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

/// RSCU result row
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RscuRow {
    pub amino_acid: String,
    pub aa_code: char,
    pub codon: String,
    pub subfam: String,
    pub count: f64,
    pub prop: f64,
    pub w_cai: f64,
    pub rscu: f64,
}

/// RSCU table
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RscuTable {
    pub rows: Vec<RscuRow>,
    /// Lookup by codon
    pub by_codon: HashMap<String, RscuRow>,
}

impl RscuTable {
    /// Get CAI weight for a codon
    pub fn w_cai(&self, codon: &str) -> f64 {
        self.by_codon
            .get(codon)
            .map(|r| r.w_cai)
            .unwrap_or(1.0)
    }

    /// Get RSCU for a codon
    pub fn rscu(&self, codon: &str) -> f64 {
        self.by_codon
            .get(codon)
            .map(|r| r.rscu)
            .unwrap_or(0.0)
    }
}

/// Estimate Relative Synonymous Codon Usage (RSCU).
///
/// # Arguments
/// * `cf` - Codon frequency matrix
/// * `weight` - Optional weights per gene (e.g., expression levels). Defaults to 1.0.
/// * `pseudo_cnt` - Pseudo count to avoid division by zero (default 1.0)
/// * `codon_table` - Genetic code table
/// * `level` - "subfam" (default) or "amino_acid"
/// * `incl_stop` - Whether to include stop codons
pub fn est_rscu(
    cf: &CodonFreqMatrix,
    weight: Option<&[f64]>,
    pseudo_cnt: f64,
    codon_table: &CodonTable,
    _level: &str,
    incl_stop: bool,
) -> RscuTable {
    let n_genes = cf.n_genes();

    // Default weight = 1.0 per gene
    let weights: Vec<f64> = if let Some(w) = weight {
        assert_eq!(w.len(), n_genes, "weights length must match number of genes");
        w.to_vec()
    } else {
        vec![1.0; n_genes]
    };

    // Compute weighted codon counts across all genes
    let mut codon_counts: HashMap<String, f64> = HashMap::new();
    for (i, row) in cf.matrix.iter().enumerate() {
        let w = weights[i];
        for (j, &count) in row.iter().enumerate() {
            if count > 0.0 {
                let codon = &cf.codons[j];
                *codon_counts.entry(codon.clone()).or_default() += count * w;
            }
        }
    }

    // Add pseudocount to all codons
    for codon in &cf.codons {
        let info = codon_table.codon_map.get(codon);
        if info.is_none() {
            continue;
        }
        let info = info.unwrap();
        if !incl_stop && info.aa_code == '*' {
            continue;
        }
        *codon_counts.entry(codon.clone()).or_default() += pseudo_cnt;
    }

    // Build groups based on analysis level
    let groups = codon_table.subfamily_groups();

    // Compute RSCU per group
    let mut rows = Vec::new();

    for (_subfam, codons_in_group) in &groups {
        // Sum of counts in this group
        let group_sum: f64 = codons_in_group
            .iter()
            .map(|c| codon_counts.get(c).copied().unwrap_or(pseudo_cnt))
            .sum();

        let n_codons = codons_in_group.len() as f64;

        for codon in codons_in_group {
            let count = codon_counts.get(codon).copied().unwrap_or(pseudo_cnt);
            let prop = if group_sum > 0.0 {
                count / group_sum
            } else {
                0.0
            };
            // RSCU = observed proportion * degeneracy = (count/sum) * n
            let rscu = if group_sum > 0.0 && n_codons > 0.0 {
                (count * n_codons) / group_sum
            } else {
                0.0
            };

            // w_cai = RSCU / max(RSCU in this group)
            let max_rscu = codons_in_group
                .iter()
                .map(|c| {
                    let ct = codon_counts.get(c).copied().unwrap_or(pseudo_cnt);
                    if group_sum > 0.0 {
                        (ct * n_codons) / group_sum
                    } else {
                        0.0
                    }
                })
                .fold(0.0f64, f64::max);

            let w_cai = if max_rscu > 0.0 { rscu / max_rscu } else { 1.0 };

            if let Some(info) = codon_table.codon_map.get(codon) {
                if !incl_stop && info.aa_code == '*' {
                    continue;
                }
                rows.push(RscuRow {
                    amino_acid: info.amino_acid.clone(),
                    aa_code: info.aa_code,
                    codon: codon.clone(),
                    subfam: info.subfam.clone(),
                    count,
                    prop,
                    w_cai,
                    rscu,
                });
            }
        }
    }

    // Sort rows by codon
    rows.sort_by(|a, b| a.codon.cmp(&b.codon));

    let by_codon: HashMap<String, RscuRow> = rows
        .iter()
        .map(|r| (r.codon.clone(), r.clone()))
        .collect();

    RscuTable { rows, by_codon }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};

    #[test]
    fn test_rscu_uniform() {
        let ct = CodonTable::standard();

        // Create a sequence with uniform codon usage
        // Use all codons equally
        let mut seq = String::new();
        for codon in &ct.all_codons {
            if !ct.is_stop(codon) {
                seq.push_str(codon);
            }
        }

        let cds = CdsSeq {
            id: "uniform".into(),
            seq: seq.as_bytes().to_vec(),
            codons: seq.as_bytes().chunks(3).map(|c| String::from_utf8_lossy(c).to_string()).collect(),
        };

        let cf = count_codons(&[cds], &ct);
        let rscu = est_rscu(&cf, None, 1.0, &ct, "subfam", false);

        // With uniform usage, all RSCU values should be close to 1.0
        for row in &rscu.rows {
            assert!((row.rscu - 1.0).abs() < 0.01,
                "RSCU for {} should be ~1.0, got {}", row.codon, row.rscu);
        }
    }

    #[test]
    fn test_rscu_with_weights() {
        let ct = CodonTable::standard();

        // Gene 1 uses only TTT (Phe)
        let cds1 = CdsSeq {
            id: "gene1".into(),
            seq: b"TTTTTT".to_vec(),
            codons: vec!["TTT".into(), "TTT".into()],
        };
        // Gene 2 uses only TTC (Phe)
        let cds2 = CdsSeq {
            id: "gene2".into(),
            seq: b"TTCTTC".to_vec(),
            codons: vec!["TTC".into(), "TTC".into()],
        };

        let cf = count_codons(&[cds1, cds2], &ct);

        // Equal weights: RSCU for TTT and TTC should both be 1.0
        let rscu_equal = est_rscu(&cf, None, 1.0, &ct, "subfam", false);
        let ttt_row = rscu_equal.by_codon.get("TTT").unwrap();
        let ttc_row = rscu_equal.by_codon.get("TTC").unwrap();
        assert!((ttt_row.rscu - 1.0).abs() < 0.01);
        assert!((ttc_row.rscu - 1.0).abs() < 0.01);

        // Weight gene1 higher: TTT should have higher RSCU
        let rscu_weighted = est_rscu(&cf, Some(&[10.0, 1.0]), 1.0, &ct, "subfam", false);
        let ttt_w = rscu_weighted.by_codon.get("TTT").unwrap();
        let ttc_w = rscu_weighted.by_codon.get("TTC").unwrap();
        assert!(ttt_w.rscu > ttc_w.rscu);
    }
}
