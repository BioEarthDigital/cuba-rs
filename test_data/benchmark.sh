#!/bin/bash
# Speed benchmark: R cubar vs Rust cuba-rs
set -e

echo "============================================"
echo "  CUBA-RS SPEED BENCHMARK: R vs Rust"
echo "============================================"
echo ""

CUBAR="target/release/cuba-rs"
DIR="test_data"

for N in 500 2000 6600; do
    FASTA="$DIR/yeast${N}.fasta"
    REF="$DIR/yeast_ref100.fasta"

    # Generate data if needed
    if [ "$N" = "2000" ] && [ ! -f "$FASTA" ]; then
        Rscript -e 'library(cubar); library(Biostrings); data(yeast_cds); writeXStringSet(yeast_cds[1:2000], "'$FASTA'", format="fasta")' 2>/dev/null
    fi
    if [ "$N" = "6600" ] && [ ! -f "$FASTA" ]; then
        Rscript -e 'library(cubar); library(Biostrings); data(yeast_cds); writeXStringSet(yeast_cds, "'$FASTA'", format="fasta")' 2>/dev/null
    fi

    [ ! -f "$FASTA" ] && FASTA="$DIR/yeast_full.fasta"

    echo "--- $N sequences ---"
    echo ""

    # R: ENC only (pure computation time)
    echo -n "R ENC computation: "
    Rscript -e '
        suppressMessages({library(cubar); library(Biostrings)})
        cds <- readDNAStringSet("'$FASTA'")
        valid <- check_cds(cds)
        cf <- count_codons(valid)
        t0 <- Sys.time()
        enc <- get_enc(cf)
        cat(sprintf("%.4f s\n", as.numeric(Sys.time()-t0, units="secs")))
    ' 2>/dev/null

    # R: Full pipeline (load + count + ENC)
    echo -n "R full pipeline:   "
    time (Rscript -e '
        suppressMessages({library(cubar); library(Biostrings)})
        cds <- readDNAStringSet("'$FASTA'")
        valid <- check_cds(cds)
        cf <- count_codons(valid)
        enc <- get_enc(cf)
        cat(sprintf("ENC range: %.2f-%.2f (mean %.2f)\n", min(enc), max(enc), mean(enc)))
    ' 2>/dev/null) 2>&1 | grep real | awk '{print $2}'

    # Rust: Full pipeline
    echo -n "Rust full pipeline: "
    time ($CUBAR enc "$FASTA" -o /dev/null 2>/dev/null) 2>&1 | grep real | awk '{print $2}'

    # R: CAI pipeline
    echo -n "R CAI:             "
    time (Rscript -e '
        suppressMessages({library(cubar); library(Biostrings)})
        cds <- readDNAStringSet("'$FASTA'")
        ref <- readDNAStringSet("'$REF'")
        valid <- check_cds(cds); ref_v <- check_cds(ref)
        cf <- count_codons(valid); ref_cf <- count_codons(ref_v)
        rscu <- est_rscu(ref_cf)
        cai <- get_cai(cf, rscu)
        cat(sprintf("CAI range: %.3f-%.3f\n", min(cai), max(cai)))
    ' 2>/dev/null) 2>&1 | grep real | awk '{print $2}'

    # Rust: CAI pipeline
    echo -n "Rust CAI:          "
    time ($CUBAR cai "$FASTA" -r "$REF" -o /dev/null 2>/dev/null) 2>&1 | grep real | awk '{print $2}'

    echo ""
done

echo "============================================"
echo "  SUMMARY TABLE"
echo "============================================"
echo ""
echo "Note: Times include I/O (FASTA reading + validation)"
echo "Times: pure computation for R, total wall-clock for both"
echo ""