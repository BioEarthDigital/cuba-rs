use std::collections::HashMap;
use rayon::prelude::*;
use crate::genetic_code::CodonTable;
use crate::sequence::CodonFreqMatrix;

/// tRNA weight data: maps anticodon -> tRNA gene copy number
pub type TrnaWeights = HashMap<String, f64>;

/// Calculate tRNA Adaptation Index (tAI) for each gene.
///
/// Reference: dos Reis M, Savva R, Wernisch L (2004) Nucleic Acids Res 32:5036-5044.
///
/// tAI measures how well a gene's codon usage matches the host's tRNA pool.
/// Values range from 0 to 1; higher = better adaptation.
///
/// # Arguments
/// * `cf` - Codon frequency matrix
/// * `trna_weights` - tRNA gene copy numbers per anticodon
/// * `codon_table` - Genetic code table
/// * `s_values` - Optional wobble pairing penalties (s_ij). If not provided, defaults are used.
pub fn get_tai(
    cf: &CodonFreqMatrix,
    trna_weights: &TrnaWeights,
    codon_table: &CodonTable,
    s_values: Option<&HashMap<String, f64>>,
) -> Vec<f64> {
    // Step 1: Compute absolute adaptiveness W_i for each codon
    let w_abs = compute_w_abs(trna_weights, s_values, codon_table);

    // Step 2: Compute relative adaptiveness w_i = W_i / max(W)
    let max_w = w_abs.values().cloned().fold(0.0f64, f64::max);

    // Compute geometric mean of non-zero w_i for fallback
    let non_zero_ws: Vec<f64> = w_abs.values().copied().filter(|&w| w > 0.0).collect();
    let geo_mean_w = if !non_zero_ws.is_empty() {
        let sum_log: f64 = non_zero_ws.iter().map(|w| w.ln()).sum();
        (sum_log / non_zero_ws.len() as f64).exp()
    } else {
        0.5
    };

    let w_rel: HashMap<String, f64> = w_abs
        .iter()
        .map(|(codon, &w)| {
            let rel = if max_w > 0.0 && w > 0.0 {
                w / max_w
            } else {
                geo_mean_w
            };
            (codon.clone(), rel)
        })
        .collect();

    // Step 3: tAI_g = geometric mean of w_i for all codons in gene
    cf.matrix
        .par_iter()
        .map(|row| {
            let mut sum_log_w = 0.0f64;
            let mut total_count = 0.0f64;
            for (j, &count) in row.iter().enumerate() {
                if count > 0.0 {
                    let codon = &cf.codons[j];
                    if codon_table.is_stop(codon) {
                        continue;
                    }
                    if let Some(&w) = w_rel.get(codon) {
                        if w > 0.0 {
                            sum_log_w += count * w.ln();
                            total_count += count;
                        }
                    }
                }
            }
            if total_count > 0.0 {
                (sum_log_w / total_count).exp()
            } else {
                0.0
            }
        })
        .collect()
}

/// Compute absolute adaptiveness W_i for each codon.
///
/// W_i = sum_j (1 - s_ij) * tGCN_ij
///
/// where tGCN_ij is the tRNA gene copy number for the j-th tRNA
/// that recognizes codon i, and s_ij is the wobble penalty.
fn compute_w_abs(
    trna_weights: &TrnaWeights,
    s_values: Option<&HashMap<String, f64>>,
    codon_table: &CodonTable,
) -> HashMap<String, f64> {
    // Map anticodon to codon based on standard wobble rules
    // For now, use a simplified mapping
    let mut w_abs: HashMap<String, f64> = HashMap::new();

    // Initialize all codons with 0
    for codon in &codon_table.all_codons {
        if !codon_table.is_stop(codon) {
            w_abs.insert(codon.clone(), 0.0);
        }
    }

    // For each anticodon with tRNA count, find which codons it recognizes
    for (anticodon, &t_gcn) in trna_weights {
        if t_gcn <= 0.0 {
            continue;
        }

        // Get the codon recognized by this anticodon
        // anticodon -> codon: complement and reverse
        let codon = anticodon_to_codon(anticodon);

        // Also consider wobble pairings
        let recognized_codons = wobble_codons(anticodon, &codon);

        for (rec_codon, wobble_type) in &recognized_codons {
            let s = if let Some(s_map) = s_values {
                s_map.get(wobble_type).copied().unwrap_or(0.5)
            } else {
                // Default s values
                match wobble_type.as_str() {
                    "WC" => 0.0,      // Watson-Crick: no penalty
                    "GU" => 0.5,      // G-U wobble
                    "IU" => 0.6,      // Inosine-U
                    "IC" => 0.6,      // Inosine-C
                    "IA" => 0.7,      // Inosine-A
                    _ => 0.5,
                }
            };

            let contrib = (1.0 - s) * t_gcn;
            *w_abs.entry(rec_codon.clone()).or_default() += contrib;
        }
    }

    w_abs
}

/// Convert a tRNA anticodon (5'->3' RNA) to the DNA codon (5'->3').
///
/// tRNA anticodon (5'→3'): CAU (Met)
/// Reverse (3'→5'): UAC
/// RNA complement → DNA: U→A, A→T, C→G = ATG (correct)
pub(crate) fn anticodon_to_codon(anticodon: &str) -> String {
    let complement = |b: char| -> char {
        match b {
            'A' => 'T',  // RNA A pairs with DNA T
            'T' => 'A',
            'U' => 'A',  // RNA U pairs with DNA A (mRNA U → tRNA A, so complement is A)
            'G' => 'C',
            'C' => 'G',
            'I' => 'A',  // Inosine pairs with A
            _ => b,
        }
    };
    anticodon.chars().rev().map(complement).collect()
}

/// Generate all codons recognized by an anticodon through wobble pairing
fn wobble_codons(anticodon: &str, primary_codon: &str) -> Vec<(String, String)> {
    let mut results = vec![(primary_codon.to_string(), "WC".to_string())];

    // The first base of the anticodon (3' end) pairs with the 3rd base of the codon
    let wobble_base = anticodon.chars().next().unwrap();

    match wobble_base {
        'G' => {
            // G can pair with C (WC) or U (wobble)
            let mut wob = primary_codon.to_string();
            wob.replace_range(2..3, "T"); // G-U wobble
            results.push((wob, "GU".to_string()));
        }
        'U' => {
            // U can pair with A (WC) or G (wobble)
            let mut wob = primary_codon.to_string();
            wob.replace_range(2..3, "G"); // U-G wobble
            results.push((wob, "GU".to_string()));
        }
        'I' => {
            // Inosine can pair with U, C, or A
            let mut wob_u = primary_codon.to_string();
            wob_u.replace_range(2..3, "T");
            results.push((wob_u, "IU".to_string()));
            let mut wob_c = primary_codon.to_string();
            wob_c.replace_range(2..3, "C");
            results.push((wob_c, "IC".to_string()));
        }
        _ => {}
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anticodon_to_codon() {
        // Anticodon 3'-UAC-5' (tRNA-Met) -> Codon 5'-ATG-3'
        let anticodon = "UAC";
        let _codon = anticodon_to_codon(anticodon);
        // U->A, A->T, C->G -> ATG... wait
        // Reverse: UAC -> CAU
        // Complement: C->G, A->T, U->A -> GTA... hmm
        // Actually:
        // Anticodon (3'->5'): U A C
        // Reverse (5'->3'): C A U
        // Complement: G T A
        // Codon: GTA
        // No wait, tRNA-Met anticodon is CAU (3'-UAC-5'), pairs with AUG (5'-AUG-3')
        // anticodon = "CAU" (5'->3')
        // complement: GTA... This doesn't seem right
        // Let me reconsider:
        // Anticodon stored in 5'->3' is CAU
        // Reverse: UAC
        // Complement: AUG -> correct!
        let anticodon2 = "CAU";
        let codon2 = anticodon_to_codon(anticodon2);
        assert_eq!(codon2, "ATG", "CAU (Met anticodon) -> ATG (Met codon)");

        // Anticodon CCA (Trp) -> Codon TGG
        // Reverse: ACC
        // Complement: TGG -> correct!
        let anticodon3 = "CCA";
        let codon3 = anticodon_to_codon(anticodon3);
        assert_eq!(codon3, "TGG", "CCA (Trp anticodon) -> TGG (Trp codon)");
    }

    #[test]
    fn test_wobble_codons() {
        // G at wobble position: pairs with C (WC) and T (wobble)
        let primary = "TTC"; // Phe
        let wobbles = wobble_codons("GAA", primary); // anticodon GAA (3'-AAG-5') for Phe
        // G pairs with C -> TTC and T -> TTT
        assert!(wobbles.iter().any(|(c, _)| c == "TTC"));
        assert!(wobbles.iter().any(|(c, _)| c == "TTT"));
    }
}
