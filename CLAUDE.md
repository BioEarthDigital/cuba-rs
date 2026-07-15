# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test

```bash
cargo build                    # debug build
cargo build --release          # release build (for benchmarking)
cargo test                     # run all tests (22 tests across cubar-core)
cargo test -p cubar-core       # run only library tests
cargo test -p cubar-core enc   # run tests matching "enc"
```

## Architecture

Cargo workspace with two crates:

**`cubar-core`** (library) — All computation logic. No CLI awareness.
- `genetic_code.rs` — `CodonTable` struct; 27 NCBI genetic codes stored as `HashMap<String, CodonInfo>`. The central data structure. `subfamily_groups()` partitions codons by first-2-bases + amino acid (e.g., Leu_TT, Leu_CT), which feeds all metrics.
- `sequence.rs` — `CdsSeq` (id + seq bytes + codon triplets), `CodonFreqMatrix` (genes × codons f64 matrix). FASTA I/O via `needletail`. `check_cds()` skips sequences with internal stop codons.
- `metrics/enc.rs` — ENC algorithm matches cubar R exactly: groups codons by subfamily, computes `p = (count+1)/(n+k)` pseudo-corrected frequencies, then `N_d = n_groups × Σn / Σ(n×f)` per degeneracy class. Sum across classes gives ENC.
- `metrics/rscu.rs` — `RscuTable` with per-codon `w_cai` (relative to max RSCU in subfamily). Supports gene-level weights and pseudocounts.
- `metrics/cai.rs` — Geometric mean of `w_cai` across all codons in a gene.
- `metrics/tai.rs` — tRNA anticodon→codon conversion (reverse + RNA→DNA complement). Wobble pairing support (G-U, Inosine).
- `optimize.rs` — `est_optimal_codons()` picks highest-count codon per subfamily. `codon_optimize()` replaces each codon with its optimal synonym.
- `slide.rs` — Sliding window over `CdsSeq.codons`, re-counts codons per window and computes metrics.

**`cubar-cli`** (binary) — Thin CLI layer. Each subcommand file in `commands/` follows the same pattern: load codon table → load FASTA → count codons → compute metric → write CSV/TSV/JSON via `write_results()`.

## Key Design Patterns

- **Data flow**: FASTA → `CdsSeq` → `CodonFreqMatrix` → metric function → output. Every command follows this pipeline.
- **Genetic code is explicit**: Every function that needs codon→AA mapping takes `&CodonTable`. No global/ambient genetic code.
- **Subfamily grouping drives everything**: `CodonTable::subfamily_groups()` returns `HashMap<String, Vec<String>>` mapping subfamily names (like "Leu_CT") to codon lists. RSCU, ENC, CAI, Fop, and optimal codon detection all iterate over these groups.
- **Stop codons are excluded per-caller**: `count_codons()` counts everything including stops; each metric's function decides whether to filter them.
- **Parity with R cubar**: The ENC algorithm was matched line-by-line to the R source. Other metrics verified against yeast_cds dataset. Any algorithm change should include a correlation check against R output.
- **rayon is used in metrics for per-gene parallelism**: `cf.matrix.par_iter().map(|row| { compute_for_gene(row) }).collect()`. Not used in GC functions (too lightweight to benefit).
