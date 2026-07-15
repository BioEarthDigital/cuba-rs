use crate::genetic_code::CodonTable;
use crate::sequence::{CdsSeq, count_codons};
use crate::metrics::{enc::get_enc, cai::get_cai, gc::get_all_gc};
use crate::metrics::rscu::RscuTable;

/// Result of a sliding window analysis for one window
#[derive(Debug, Clone, serde::Serialize)]
pub struct SlideWindow {
    pub seq_id: String,
    pub start: usize,      // 1-indexed start position
    pub end: usize,        // 1-indexed end position
    pub window_index: usize,
    pub enc: Option<f64>,
    pub cai: Option<f64>,
    pub gc: Option<f64>,
    pub gc3s: Option<f64>,
    pub gc4d: Option<f64>,
}

/// Run a sliding window analysis on a set of coding sequences.
///
/// # Arguments
/// * `cds_seqs` - Coding sequences
/// * `window_size` - Window size in codons
/// * `step_size` - Step size in codons
/// * `codon_table` - Genetic code table
/// * `rscu` - Optional RSCU table for CAI calculation
/// * `metrics` - Which metrics to compute ("enc", "cai", "gc", or "all")
pub fn slide(
    cds_seqs: &[CdsSeq],
    window_size: usize,
    step_size: usize,
    codon_table: &CodonTable,
    rscu: Option<&RscuTable>,
    metrics: &[&str],
) -> Vec<SlideWindow> {
    let mut results = Vec::new();
    let compute_enc = metrics.contains(&"enc") || metrics.contains(&"all");
    let compute_cai = metrics.contains(&"cai") || metrics.contains(&"all");
    let compute_gc = metrics.contains(&"gc") || metrics.contains(&"all");

    for seq in cds_seqs {
        let n_codons = seq.codons.len();
        if n_codons < window_size {
            continue;
        }

        let mut window_idx = 0;
        let mut pos = 0;

        while pos + window_size <= n_codons {
            let window_codons = &seq.codons[pos..pos + window_size];

            // Create a temporary CdsSeq for this window
            let window_seq = CdsSeq {
                id: format!("{}_w{}", seq.id, window_idx),
                seq: window_codons.join("").as_bytes().to_vec(),
                codons: window_codons.to_vec(),
            };

            let cf = count_codons(&[window_seq], codon_table);

            let enc = if compute_enc {
                Some(get_enc(&cf, codon_table, "subfam")[0])
            } else {
                None
            };

            let cai = if compute_cai {
                rscu.map(|r| get_cai(&cf, r, "subfam")[0])
            } else {
                None
            };

            let (gc, gc3s, gc4d) = if compute_gc {
                let gc_results = &get_all_gc(&cf, codon_table)[0];
                (Some(gc_results.gc), Some(gc_results.gc3s), Some(gc_results.gc4d))
            } else {
                (None, None, None)
            };

            results.push(SlideWindow {
                seq_id: seq.id.clone(),
                start: pos * 3 + 1,
                end: (pos + window_size) * 3,
                window_index: window_idx,
                enc,
                cai,
                gc,
                gc3s,
                gc4d,
            });

            window_idx += 1;
            pos += step_size;
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic_code::CodonTable;

    #[test]
    fn test_slide() {
        let ct = CodonTable::standard();

        // Create a long CDS with 10 codons
        let codons: Vec<String> = vec![
            "ATG", "GCT", "GCC", "GCA", "GCG",
            "TGG", "GGT", "GGC", "GGA", "TAA",
        ].iter().map(|s| s.to_string()).collect();

        let seq_str = codons.join("");
        let cds = CdsSeq {
            id: "test".into(),
            seq: seq_str.as_bytes().to_vec(),
            codons,
        };

        let results = slide(&[cds], 5, 3, &ct, None, &["enc", "gc"]);

        // With window=5, step=3, and 10 codons:
        // Window 0: codons 0-4, Window 1: codons 3-7, Window 2: codons 6-9 (no, 6+5=11>10, so only 2)
        // Actually: 10 codons, window=5, step=3
        // pos=0: 0..5 ✓, pos=3: 3..8 ✓, pos=6: 6..11 ✗ (no: 10 codons, pos+window=6+5=11>10)
        // Wait: pos=3: 3..8 ✓ (8<10), pos=6: 6..11 ✗
        // Oh wait, we have 10 codons (indices 0-9). pos+window <= n_codons:
        // pos=0: 0+5=5 <= 10 ✓
        // pos=3: 3+5=8 <= 10 ✓
        // pos=6: 6+5=11 <= 10 ✗
        // So 2 windows
        assert_eq!(results.len(), 2);

        // Check first window
        assert_eq!(results[0].seq_id, "test");
        assert_eq!(results[0].start, 1);
        assert_eq!(results[0].end, 15);
        assert!(results[0].enc.is_some());
        assert!(results[0].gc.is_some());

        // Check second window
        assert_eq!(results[1].start, 10);
        assert_eq!(results[1].end, 24);
    }
}
