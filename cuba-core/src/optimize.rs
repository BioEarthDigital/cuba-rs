use std::collections::HashMap;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;
use crate::metrics::rscu::est_rscu;

/// Table of optimal codons per amino acid
#[derive(Debug, Clone)]
pub struct OptimalCodonTable {
    /// amino acid 1-letter -> optimal codon
    pub optimal: HashMap<char, String>,
    /// amino acid 1-letter -> all optimal codons with their preferences
    pub optimal_list: HashMap<char, Vec<String>>,
}

/// Estimate optimal codons based on codon frequencies.
///
/// Optimal codons are those that:
/// 1. Increase in frequency with increasing gene expression
/// 2. Have highest RSCU in highly expressed genes
pub fn est_optimal_codons(
    cf: &CodonFreqMatrix,
    codon_table: &CodonTable,
) -> OptimalCodonTable {
    // Simple approach: for each subfamily, the codon with highest overall usage
    let rscu = est_rscu(cf, None, 1.0, codon_table, "subfam", false);
    let groups = codon_table.subfamily_groups();

    let mut optimal: HashMap<char, String> = HashMap::new();
    let mut optimal_list: HashMap<char, Vec<String>> = HashMap::new();

    for (_sf, codons) in &groups {
        if codons.is_empty() {
            continue;
        }

        let aa = codon_table.codon_map.get(&codons[0]).map(|ci| ci.aa_code);
        if let Some(aa) = aa {
            if aa == '*' {
                continue;
            }

            // Find codon with highest count in this subfamily
            let best_codon = codons
                .iter()
                .max_by(|a, b| {
                    let count_a = rscu.by_codon.get(*a).map(|r| r.count).unwrap_or(0.0);
                    let count_b = rscu.by_codon.get(*b).map(|r| r.count).unwrap_or(0.0);
                    count_a.partial_cmp(&count_b).unwrap_or(std::cmp::Ordering::Equal)
                })
                .cloned();

            if let Some(best) = best_codon {
                // For the "optimal" table, use the single best codon per amino acid
                // (or the best for each subfamily)
                optimal.entry(aa).or_insert_with(|| best.clone());
                optimal_list.entry(aa).or_default().push(best);
            }
        }
    }

    OptimalCodonTable {
        optimal,
        optimal_list,
    }
}

/// Optimize codon usage of a DNA sequence.
///
/// Replaces each codon with the corresponding synonymous optimal codon.
///
/// # Arguments
/// * `seq` - DNA sequence as string
/// * `optimal_codons` - Table of optimal codons
/// * `method` - Currently only "naive" is supported
pub fn codon_optimize(
    seq: &str,
    optimal_codons: &OptimalCodonTable,
    _method: &str,
) -> String {
    let seq = seq.to_uppercase();
    let mut result = String::with_capacity(seq.len());
    let mut i = 0;
    let bytes = seq.as_bytes();

    while i + 3 <= bytes.len() {
        let codon = &seq[i..i + 3];

        // Find the amino acid for this codon
        // (We need a codon table lookup)
        // For now, use the optimal table directly
        let optimized = if let Some(aa) = codon_to_aa(codon) {
            if let Some(opt_codon) = optimal_codons.optimal.get(&aa) {
                opt_codon.as_str()
            } else {
                codon
            }
        } else {
            codon
        };

        result.push_str(optimized);
        i += 3;
    }

    // Append any remaining bases
    if i < bytes.len() {
        result.push_str(&seq[i..]);
    }

    result
}

/// Quick codon-to-amino-acid lookup (standard code)
fn codon_to_aa(codon: &str) -> Option<char> {
    match codon {
        "TTT" | "TTC" => Some('F'),
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => Some('L'),
        "ATT" | "ATC" | "ATA" => Some('I'),
        "ATG" => Some('M'),
        "GTT" | "GTC" | "GTA" | "GTG" => Some('V'),
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => Some('S'),
        "CCT" | "CCC" | "CCA" | "CCG" => Some('P'),
        "ACT" | "ACC" | "ACA" | "ACG" => Some('T'),
        "GCT" | "GCC" | "GCA" | "GCG" => Some('A'),
        "TAT" | "TAC" => Some('Y'),
        "CAT" | "CAC" => Some('H'),
        "CAA" | "CAG" => Some('Q'),
        "AAT" | "AAC" => Some('N'),
        "AAA" | "AAG" => Some('K'),
        "GAT" | "GAC" => Some('D'),
        "GAA" | "GAG" => Some('E'),
        "TGT" | "TGC" => Some('C'),
        "TGG" => Some('W'),
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Some('R'),
        "GGT" | "GGC" | "GGA" | "GGG" => Some('G'),
        "TAA" | "TAG" | "TGA" => Some('*'),
        _ => None,
    }
}

/// Identify optimal codons using gene expression data.
///
/// Codons whose frequency increases with expression level are considered optimal.
pub fn est_optimal_codons_with_expression(
    cf: &CodonFreqMatrix,
    expression: &[f64],
    codon_table: &CodonTable,
) -> OptimalCodonTable {
    // Weight RSCU by expression level
    let rscu = est_rscu(cf, Some(expression), 1.0, codon_table, "subfam", false);
    let groups = codon_table.subfamily_groups();

    let mut optimal: HashMap<char, String> = HashMap::new();
    let mut optimal_list: HashMap<char, Vec<String>> = HashMap::new();

    for (_sf, codons) in &groups {
        if codons.is_empty() {
            continue;
        }

        let aa = codon_table.codon_map.get(&codons[0]).map(|ci| ci.aa_code);
        if let Some(aa) = aa {
            if aa == '*' {
                continue;
            }

            // Find codon with highest RSCU (weighted by expression)
            let best_codon = codons
                .iter()
                .max_by(|a, b| {
                    let rscu_a = rscu.by_codon.get(*a).map(|r| r.rscu).unwrap_or(0.0);
                    let rscu_b = rscu.by_codon.get(*b).map(|r| r.rscu).unwrap_or(0.0);
                    rscu_a.partial_cmp(&rscu_b).unwrap_or(std::cmp::Ordering::Equal)
                })
                .cloned();

            if let Some(best) = best_codon {
                optimal.entry(aa).or_insert_with(|| best.clone());
                optimal_list.entry(aa).or_default().push(best);
            }
        }
    }

    OptimalCodonTable {
        optimal,
        optimal_list,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;
    use crate::sequence::{CdsSeq, count_codons};

    #[test]
    fn test_codon_optimize() {
        let _ct = CodonTable::standard();

        // Create optimal codon table: prefer GCT for Ala
        let mut optimal = HashMap::new();
        optimal.insert('A', "GCT".to_string());
        let mut optimal_list = HashMap::new();
        optimal_list.insert('A', vec!["GCT".to_string()]);
        let opt_table = OptimalCodonTable { optimal, optimal_list };

        // GCC (Ala) -> GCT (optimal Ala)
        let seq = "GCC";
        let optimized = codon_optimize(seq, &opt_table, "naive");
        assert_eq!(optimized, "GCT");
    }

    #[test]
    fn test_est_optimal_codons() {
        let ct = CodonTable::standard();

        // Gene using only GCT for Ala
        let cds = CdsSeq {
            id: "test".into(),
            seq: b"GCTGCTGCT".to_vec(),
            codons: vec!["GCT".into(), "GCT".into(), "GCT".into()],
        };

        let cf = count_codons(&[cds], &ct);
        let optimal = est_optimal_codons(&cf, &ct);

        // GCT should be optimal for Ala
        assert_eq!(optimal.optimal.get(&'A'), Some(&"GCT".to_string()));
    }
}
