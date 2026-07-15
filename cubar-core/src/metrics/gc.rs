use std::collections::HashMap;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

/// GC content metrics for a gene
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GcResults {
    pub gene_id: String,
    /// Overall GC content
    pub gc: f64,
    /// GC content at the 3rd synonymous codon positions
    pub gc3s: f64,
    /// GC content at the 4-fold degenerate 3rd positions
    pub gc4d: f64,
}

/// Calculate overall GC content for each gene.
pub fn get_gc(cf: &CodonFreqMatrix) -> Vec<f64> {
    // GC content = (G + C) / total
    cf.matrix
        .iter()
        .map(|row| {
            let mut gc = 0.0f64;
            let mut total = 0.0f64;
            for (j, &count) in row.iter().enumerate() {
                if count > 0.0 {
                    let codon = &cf.codons[j];
                    let gc_count = codon.chars().filter(|&c| c == 'G' || c == 'C').count() as f64;
                    gc += gc_count * count;
                    total += count * 3.0; // 3 bases per codon
                }
            }
            if total > 0.0 { gc / total } else { 0.0 }
        })
        .collect()
}

/// Calculate GC3s: GC content at 3rd synonymous codon positions.
///
/// 3rd synonymous positions are positions where any nucleotide change results
/// in a synonymous substitution (same amino acid).
///
/// In practice, this includes codons from 4-fold degenerate families
/// (Ala, Gly, Pro, Thr, Val) and some others.
pub fn get_gc3s(cf: &CodonFreqMatrix, codon_table: &CodonTable) -> Vec<f64> {
    let groups = codon_table.subfamily_groups();

    // Identify synonymous 3rd positions
    // These are positions within subfamilies where synonymous codons differ only at position 3
    let mut syn_positions: HashMap<String, bool> = HashMap::new();

    for (_sf, codons) in &groups {
        if codons.len() >= 2 {
            // Check if codons in this subfamily differ only in position 3
            let first_two: Vec<&str> = codons.iter().map(|c| &c[..2]).collect();
            let all_same_prefix = first_two.windows(2).all(|w| w[0] == w[1]);
            if all_same_prefix {
                // All codons in this family differ only at position 3
                // The 3rd position is synonymous
                for codon in codons {
                    if !codon_table.is_stop(codon) {
                        syn_positions.insert(codon.clone(), true);
                    }
                }
            }
        }
    }

    cf.matrix
        .iter()
        .map(|row| {
            let mut gc = 0.0f64;
            let mut total = 0.0f64;
            for (j, &count) in row.iter().enumerate() {
                if count > 0.0 {
                    let codon = &cf.codons[j];
                    if syn_positions.contains_key(codon) {
                        let third_base = codon.chars().nth(2).unwrap();
                        if third_base == 'G' || third_base == 'C' {
                            gc += count;
                        }
                        total += count;
                    }
                }
            }
            if total > 0.0 { gc / total } else { 0.0 }
        })
        .collect()
}

/// Calculate GC4d: GC content at 4-fold degenerate 3rd positions.
///
/// 4-fold degenerate positions are those where any of the 4 nucleotides
/// at the 3rd position encodes the same amino acid.
/// This includes: Ala (GCN), Gly (GGN), Pro (CCN), Thr (ACN), Val (GTN),
/// and the 4-fold parts of Arg (CGN), Leu (CTN).
pub fn get_gc4d(cf: &CodonFreqMatrix, codon_table: &CodonTable) -> Vec<f64> {
    let groups = codon_table.subfamily_groups();

    // Identify 4-fold degenerate subfamilies: exactly 4 codons, same first 2 bases
    let mut fourfold_codons: HashMap<String, bool> = HashMap::new();

    for (_sf, codons) in &groups {
        if codons.len() == 4 {
            let first_two: Vec<&str> = codons.iter().map(|c| &c[..2]).collect();
            if first_two.windows(2).all(|w| w[0] == w[1]) {
                for codon in codons {
                    if !codon_table.is_stop(codon) {
                        fourfold_codons.insert(codon.clone(), true);
                    }
                }
            }
        }
    }

    cf.matrix
        .iter()
        .map(|row| {
            let mut gc = 0.0f64;
            let mut total = 0.0f64;
            for (j, &count) in row.iter().enumerate() {
                if count > 0.0 {
                    let codon = &cf.codons[j];
                    if fourfold_codons.contains_key(codon) {
                        let third_base = codon.chars().nth(2).unwrap();
                        if third_base == 'G' || third_base == 'C' {
                            gc += count;
                        }
                        total += count;
                    }
                }
            }
            if total > 0.0 { gc / total } else { 0.0 }
        })
        .collect()
}

/// Get all GC metrics at once.
pub fn get_all_gc(cf: &CodonFreqMatrix, codon_table: &CodonTable) -> Vec<GcResults> {
    let gc = get_gc(cf);
    let gc3s = get_gc3s(cf, codon_table);
    let gc4d = get_gc4d(cf, codon_table);

    cf.gene_ids
        .iter()
        .enumerate()
        .map(|(i, id)| GcResults {
            gene_id: id.clone(),
            gc: gc[i],
            gc3s: gc3s[i],
            gc4d: gc4d[i],
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};

    #[test]
    fn test_gc() {
        let ct = CodonTable::standard();
        // ATG (A,T,G) + CCC (C,C,C) = 6 bases, G=1, C=3, GC=4/6=0.6667
        let cds = CdsSeq {
            id: "test".into(),
            seq: b"ATGCCC".to_vec(),
            codons: vec!["ATG".into(), "CCC".into()],
        };
        let cf = count_codons(&[cds], &ct);
        let gc_vals = get_gc(&cf);
        assert!((gc_vals[0] - 0.6667).abs() < 0.01, "GC = {:.4}", gc_vals[0]);
    }

    #[test]
    fn test_gc4d() {
        let ct = CodonTable::standard();
        // GCC -> Ala, 4-fold degenerate, 3rd base is C (GC)
        // GCA -> Ala, 4-fold degenerate, 3rd base is A (not GC)
        let cds = CdsSeq {
            id: "test".into(),
            seq: b"GCCGCA".to_vec(),
            codons: vec!["GCC".into(), "GCA".into()],
        };
        let cf = count_codons(&[cds], &ct);
        let gc4d_vals = get_gc4d(&cf, &ct);
        assert!((gc4d_vals[0] - 0.5).abs() < 0.01, "GC4d = {:.4}", gc4d_vals[0]);
    }
}
