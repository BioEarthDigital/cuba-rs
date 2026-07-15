use rayon::prelude::*;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;
use crate::metrics::rscu::RscuTable;

/// Calculate Codon Adaptation Index (CAI) for each gene.
///
/// CAI measures how well a gene's codon usage matches that of highly expressed genes.
/// Values range from 0 to 1. Higher values indicate better adaptation.
///
/// Reference: Sharp PM, Li WH (1987) Nucleic Acids Res 15:1281-1295.
///
/// # Arguments
/// * `cf` - Codon frequency matrix for target genes
/// * `rscu` - RSCU table from reference genes (highly expressed), containing `w_cai` weights
/// * `level` - "subfam" (default) or "amino_acid"
///
/// # Returns
/// Vector of CAI values, one per gene.
pub fn get_cai(cf: &CodonFreqMatrix, rscu: &RscuTable, level: &str) -> Vec<f64> {
    // CAI = geometric mean of w_cai for each codon in the gene
    cf.matrix
        .par_iter()
        .map(|row| {
            compute_cai_gene(row, &cf.codons, rscu, level)
        })
        .collect()
}

fn compute_cai_gene(
    row: &[f64],
    codon_names: &[String],
    rscu: &RscuTable,
    _level: &str,
) -> f64 {
    // Collect the w_cai values for each codon occurrence
    // CAI = exp( (1/L) * sum(ln(w_i)) )
    let mut sum_log_w = 0.0f64;
    let mut total_count = 0.0f64;

    for (j, &count) in row.iter().enumerate() {
        if count > 0.0 {
            let codon = &codon_names[j];
            // Skip Met, Trp, and stop codons (non-degenerate or stop)
            if let Some(rscu_row) = rscu.by_codon.get(codon) {
                let w = rscu_row.w_cai;
                if w > 0.0 {
                    sum_log_w += count * w.ln();
                    total_count += count;
                }
            }
        }
    }

    if total_count == 0.0 {
        return 0.0;
    }

    (sum_log_w / total_count).exp()
}

/// Calculate CAI using reference codon frequencies directly (alternative method).
///
/// Uses reference codon frequencies to compute w_cai internally.
pub fn get_cai_from_ref(
    cf: &CodonFreqMatrix,
    ref_cf: &CodonFreqMatrix,
    codon_table: &CodonTable,
    level: &str,
) -> Vec<f64> {
    // Compute RSCU from reference genes
    let rscu = crate::metrics::rscu::est_rscu(ref_cf, None, 1.0, codon_table, level, false);
    get_cai(cf, &rscu, level)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};
    use crate::metrics::rscu::est_rscu;

    #[test]
    fn test_cai_identical_ref() {
        let ct = CodonTable::standard();

        // Use same gene as reference and target
        let groups = ct.subfamily_groups();
        let mut seq = String::new();
        for (_sf, codons) in &groups {
            for codon in codons {
                seq.push_str(codon);
            }
        }

        let cds = CdsSeq {
            id: "test".into(),
            seq: seq.as_bytes().to_vec(),
            codons: seq.as_bytes().chunks(3).map(|c| String::from_utf8_lossy(c).to_string()).collect(),
        };

        let cf = count_codons(&[cds.clone()], &ct);
        let rscu = est_rscu(&cf, None, 1.0, &ct, "subfam", false);
        let cai = get_cai(&cf, &rscu, "subfam");

        // When target == reference, CAI should be close to 1.0 for each codon group
        assert!(cai[0] > 0.9, "CAI should be high when gene matches reference, got {}", cai[0]);
    }

    #[test]
    fn test_cai_range() {
        let ct = CodonTable::standard();

        // Reference: only uses GCT for Ala
        let ref_cds = CdsSeq {
            id: "ref".into(),
            seq: b"GCTGCTGCT".to_vec(),
            codons: vec!["GCT".into(), "GCT".into(), "GCT".into()],
        };

        // Target: uses GCT for Ala (same as ref)
        let target_cds = CdsSeq {
            id: "target".into(),
            seq: b"GCTGCTGCT".to_vec(),
            codons: vec!["GCT".into(), "GCT".into(), "GCT".into()],
        };

        let ref_cf = count_codons(&[ref_cds], &ct);
        let target_cf = count_codons(&[target_cds], &ct);
        let rscu = est_rscu(&ref_cf, None, 1.0, &ct, "subfam", false);
        let cai = get_cai(&target_cf, &rscu, "subfam");

        assert!(cai[0] <= 1.0, "CAI should be <= 1.0");
        assert!(cai[0] >= 0.0, "CAI should be >= 0.0");
    }
}
