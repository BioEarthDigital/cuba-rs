# R comparison script — computes metrics using cubar package
library(cubar)
library(Biostrings)

args <- commandArgs(trailingOnly=TRUE)
fasta_file <- if (length(args) >= 1) args[1] else "test_data/yeast500.fasta"
out_dir <- if (length(args) >= 2) args[2] else "test_data"

cat("Loading CDS from", fasta_file, "\n")
cds <- readDNAStringSet(fasta_file)
cat("Sequences:", length(cds), "\n")

# Validate CDS
valid_cds <- check_cds(cds)
cat("Valid CDS:", length(valid_cds), "\n")

# Count codons
cat("\n--- Counting codons ---\n")
t0 <- Sys.time()
cf <- count_codons(valid_cds)
t1 <- Sys.time()
cat(sprintf("Time: %.3f s\n", as.numeric(t1 - t0, units="secs")))

# ENC
cat("\n--- ENC ---\n")
t0 <- Sys.time()
enc <- get_enc(cf)
t1 <- Sys.time()
cat(sprintf("Time: %.3f s\n", as.numeric(t1 - t0, units="secs")))
cat(sprintf("ENC range: %.3f - %.3f, mean: %.3f, median: %.3f\n",
    min(enc), max(enc), mean(enc), median(enc)))

# RSCU
cat("\n--- RSCU ---\n")
t0 <- Sys.time()
rscu <- est_rscu(cf)
t1 <- Sys.time()
cat(sprintf("Time: %.3f s\n", as.numeric(t1 - t0, units="secs")))

# CAI (use top 100 highly expressed as reference — simulate from first 100 genes)
cat("\n--- CAI ---\n")
ref_genes <- head(rownames(cf), 100)
ref_cf <- cf[ref_genes, , drop=FALSE]
rscu_ref <- est_rscu(ref_cf)
t0 <- Sys.time()
cai <- get_cai(cf, rscu_ref)
t1 <- Sys.time()
cat(sprintf("Time: %.3f s\n", as.numeric(t1 - t0, units="secs")))
cat(sprintf("CAI range: %.3f - %.3f, mean: %.3f, median: %.3f\n",
    min(cai), max(cai), mean(cai), median(cai)))

# GC metrics
cat("\n--- GC ---\n")
t0 <- Sys.time()
gc_val <- get_gc(cf)
gc3s_val <- get_gc3s(cf)
gc4d_val <- get_gc4d(cf)
t1 <- Sys.time()
cat(sprintf("Time: %.3f s\n", as.numeric(t1 - t0, units="secs")))
cat(sprintf("GC range: %.3f - %.3f\n", min(gc_val), max(gc_val)))
cat(sprintf("GC3s range: %.3f - %.3f\n", min(gc3s_val), max(gc3s_val)))
cat(sprintf("GC4d range: %.3f - %.3f\n", min(gc4d_val), max(gc4d_val)))

# Save results for comparison
write.csv(data.frame(gene_id=names(enc), enc=enc),
    file.path(out_dir, "r_enc.csv"), row.names=FALSE, quote=FALSE)
write.csv(data.frame(gene_id=names(cai), cai=cai),
    file.path(out_dir, "r_cai.csv"), row.names=FALSE, quote=FALSE)
write.csv(data.frame(gene_id=names(gc_val), gc=gc_val, gc3s=gc3s_val, gc4d=gc4d_val),
    file.path(out_dir, "r_gc.csv"), row.names=FALSE, quote=FALSE)
write.csv(as.data.frame(rscu),
    file.path(out_dir, "r_rscu.csv"), row.names=FALSE, quote=FALSE)

# Write codon count summary
gene_totals <- rowSums(cf)
write.csv(data.frame(gene_id=names(gene_totals), n_codons=as.vector(gene_totals)),
    file.path(out_dir, "r_codon_counts.csv"), row.names=FALSE, quote=FALSE)

cat("\nR results saved to", out_dir, "\n")
