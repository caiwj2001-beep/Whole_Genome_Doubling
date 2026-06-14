# ============================================================
# Phase 2: Genomic Correlates of WGD 
# TP53, RB1 mutations, TMB, MSI, Mutational Signatures
# ============================================================

.libPaths(c("/home/caiwj2001/R/library", .libPaths()))
suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

DATA_DIR <- "/home/caiwj2001/肿瘤突变特征与治疗反应/msk_impact_50k_2026"
OUT_DIR  <- "/home/caiwj2001/癌症演化"

# Load merged sample data
sam <- readRDS(file.path(OUT_DIR, "data", "sample_merged.rds"))
cat("Loaded", nrow(sam), "samples\n")

# ---- 1. Load mutation data ----
cat("\n=== Loading mutation data ===\n")
mut <- fread(file.path(DATA_DIR, "data_mutations.txt"))
cat("Mutations:", nrow(mut), "\n")
cat("Columns:", paste(names(mut)[1:10], collapse=", "), "...\n")

# Identify key columns
cat("Sample ID column check:\n")
print(head(mut$Tumor_Sample_Barcode, 3))

# ---- 2. Build gene-level binary mutation matrix ----
# Focus on key cancer genes
KEY_GENES <- c("TP53", "RB1", "KRAS", "EGFR", "PIK3CA", "PTEN", "APC", "BRAF", 
               "CDKN2A", "NF1", "ARID1A", "CTNNB1", "SMAD4", "ATM", "BRCA1", 
               "BRCA2", "KMT2D", "KMT2C", "FAT1", "NOTCH1", "ERBB2", "IDH1",
               "KEAP1", "STK11", "NFE2L2", "RICTOR", "TSC1", "TSC2")

# Get sample IDs from mutations
mut$Tumor_Sample_Barcode <- as.character(mut$Tumor_Sample_Barcode)

# Binary mutation per gene
gene_mut_list <- list()
for (gene in KEY_GENES) {
  samples_with_mut <- unique(mut[Hugo_Symbol == gene, Tumor_Sample_Barcode])
  gene_mut_list[[gene]] <- samples_with_mut
}

# Build matrix
all_samples <- unique(sam$SAMPLE_ID)
mut_mat <- data.frame(SAMPLE_ID = all_samples, stringsAsFactors = FALSE)
for (gene in KEY_GENES) {
  mut_mat[[paste0("MUT_", gene)]] <- ifelse(mut_mat$SAMPLE_ID %in% gene_mut_list[[gene]], 1, 0)
}

cat("\nMutation matrix built:", nrow(mut_mat), "samples x", length(KEY_GENES), "genes\n")
cat("TP53 mutated:", sum(mut_mat$MUT_TP53), "\n")
cat("RB1 mutated:", sum(mut_mat$MUT_RB1), "\n")

# ---- 3. Merge with sample data ----
sam <- merge(sam, mut_mat, by = "SAMPLE_ID", all.x = TRUE)
# Fill NAs with 0 for samples not in mutation data
for (gene in KEY_GENES) {
  col_name <- paste0("MUT_", gene)
  sam[[col_name]][is.na(sam[[col_name]])] <- 0
}

# ---- 4. Association: WGD vs Key Driver Mutations ----
cat("\n=== WGD vs Driver Mutations (Fisher's exact test) ===\n")
assoc_results <- data.frame(Gene = character(), OR = numeric(), 
                            CI_Lower = numeric(), CI_Upper = numeric(),
                            P_value = numeric(), WGDplus_MutRate = numeric(),
                            WGDminus_MutRate = numeric(), stringsAsFactors = FALSE)

for (gene in KEY_GENES) {
  col_name <- paste0("MUT_", gene)
  tbl <- table(sam$WGD, sam[[col_name]])
  if (nrow(tbl) >= 2 && ncol(tbl) >= 2) {
    ft <- fisher.test(tbl)
    wgd_plus_rate <- 100 * tbl[2,2] / sum(tbl[2,])
    wgd_minus_rate <- 100 * tbl[1,2] / sum(tbl[1,])
    assoc_results <- rbind(assoc_results, data.frame(
      Gene = gene, OR = round(ft$estimate, 2),
      CI_Lower = round(ft$conf.int[1], 2),
      CI_Upper = round(ft$conf.int[2], 2),
      P_value = ft$p.value,
      WGDplus_MutRate = round(wgd_plus_rate, 1),
      WGDminus_MutRate = round(wgd_minus_rate, 1),
      stringsAsFactors = FALSE
    ))
  }
}

# FDR correction
assoc_results$FDR <- p.adjust(assoc_results$P_value, method = "BH")
assoc_results <- assoc_results[order(assoc_results$FDR), ]
print(assoc_results, nrows = 30, na.print = "")

# ---- 5. Save association results ----
write.csv(assoc_results, file.path(OUT_DIR, "tables", "TableS2_WGD_DriverGene_Associations.csv"), row.names = FALSE)

# ---- 6. Combined TP53/RB1 analysis ----
sam$TP53_RB1_status <- "Neither"
sam$TP53_RB1_status[sam$MUT_TP53 == 1 & sam$MUT_RB1 == 0] <- "TP53_only"
sam$TP53_RB1_status[sam$MUT_TP53 == 0 & sam$MUT_RB1 == 1] <- "RB1_only"
sam$TP53_RB1_status[sam$MUT_TP53 == 1 & sam$MUT_RB1 == 1] <- "Both"
sam$TP53_RB1_status <- factor(sam$TP53_RB1_status, 
                               levels = c("Neither", "TP53_only", "RB1_only", "Both"))

cat("\n=== TP53/RB1 combined WGD rates ===\n")
tp53rb1 <- sam[, .(N = .N, WGD_Rate = round(100 * mean(WGD, na.rm = TRUE), 1)),
               by = TP53_RB1_status]
print(tp53rb1)

# ---- 7. Load and merge mutational signatures ----
cat("\n=== Loading mutational signatures ===\n")
sig_raw <- fread(file.path(DATA_DIR, "data_mutational_signatures_contribution_v2.txt"), header = TRUE, sep = "\t")
cat("Signature rows:", nrow(sig_raw), "\n")
# Print first few entity IDs
print(head(sig_raw$ENTITY_STABLE_ID, 10))

# Build sample x signature matrix
sig_ids <- sig_raw$ENTITY_STABLE_ID
sample_cols <- setdiff(names(sig_raw), c("ENTITY_STABLE_ID", "NAME", "DESCRIPTION", "URL"))
sig_mat <- as.data.frame(sig_raw[, ..sample_cols])
rownames(sig_mat) <- sig_ids
sig_df <- as.data.frame(t(sig_mat))
sig_df$SAMPLE_ID <- rownames(sig_df)

# Rename columns to remove prefix
colnames(sig_df) <- gsub("mutational_signatures_contribution_v2_", "", colnames(sig_df))
names(sig_df)[names(sig_df) == "ENTITY_STABLE_ID"] <- "SIG_ID"

cat("Signature matrix:", nrow(sig_df), "samples x", ncol(sig_df)-1, "signatures\n")

# Merge with sample data
sam <- merge(sam, sig_df, by = "SAMPLE_ID", all.x = TRUE)
cat("After signature merge:", nrow(sam), "samples\n")

# ---- 8. Signature groups ----
# COSMIC v2 signatures
SIG_APOBEC <- c("Signature2", "Signature13")
SIG_HRD    <- c("Signature3")
SIG_MMR    <- c("Signature6", "Signature15", "Signature20", "Signature26")
SIG_POLE   <- c("Signature10")
SIG_SMOKING <- c("Signature4")
SIG_AGING  <- c("Signature1")
SIG_UV     <- c("Signature7")

# Compute composite scores
available_sigs <- intersect(names(sam), c(SIG_APOBEC, SIG_HRD, SIG_MMR, SIG_POLE, SIG_SMOKING, SIG_AGING, SIG_UV))
cat("Available signatures:", paste(available_sigs, collapse=", "), "\n")

# Convert to numeric
for (s in available_sigs) {
  sam[[s]] <- as.numeric(as.character(sam[[s]]))
}

# Sum scores
APOBEC_sigs <- intersect(SIG_APOBEC, names(sam))
MMR_sigs <- intersect(SIG_MMR, names(sam))

if (length(APOBEC_sigs) > 0) sam$SIG_APOBEC <- rowSums(sam[, ..APOBEC_sigs], na.rm = TRUE)
if (length(MMR_sigs) > 0) sam$SIG_MMR <- rowSums(sam[, ..MMR_sigs], na.rm = TRUE)
if ("Signature3" %in% names(sam)) sam$SIG_HRD <- as.numeric(sam$Signature3)
if ("Signature1" %in% names(sam)) sam$SIG_AGING <- as.numeric(sam$Signature1)
if ("Signature4" %in% names(sam)) sam$SIG_SMOKING <- as.numeric(sam$Signature4)

# ---- 9. WGD vs mutational signatures ----
cat("\n=== WGD vs Mutational Signatures (Wilcoxon) ===\n")
sig_cols <- c("SIG_APOBEC", "SIG_HRD", "SIG_MMR", "SIG_AGING", "SIG_SMOKING")
sig_cols <- intersect(sig_cols, names(sam))

sig_assoc <- data.frame(Signature = character(), WGDplus_Median = numeric(),
                        WGDminus_Median = numeric(), P_value = numeric(),
                        stringsAsFactors = FALSE)

for (sc in sig_cols) {
  wgd_p <- sam[WGD == 1 & !is.na(get(sc)), median(get(sc), na.rm = TRUE)]
  wgd_m <- sam[WGD == 0 & !is.na(get(sc)), median(get(sc), na.rm = TRUE)]
  wt <- wilcox.test(sam[WGD == 1, get(sc)], sam[WGD == 0, get(sc)])
  sig_assoc <- rbind(sig_assoc, data.frame(
    Signature = sc, WGDplus_Median = round(wgd_p, 2),
    WGDminus_Median = round(wgd_m, 2), P_value = wt$p.value,
    stringsAsFactors = FALSE
  ))
}

sig_assoc$FDR <- p.adjust(sig_assoc$P_value, method = "BH")
print(sig_assoc)

# ---- 10. Save expanded dataset ----
saveRDS(sam, file.path(OUT_DIR, "data", "sample_with_genomics.rds"))
cat("\nSaved expanded dataset\n")

# ---- 11. Key summary statistics for manuscript ----
cat("\n========================================\n")
cat("KEY STATISTICS FOR MANUSCRIPT\n")
cat("========================================\n")
cat("Total tumors with WGD data:", sum(!is.na(sam$WGD)), "\n")
cat("WGD+ tumors:", sum(sam$WGD == 1, na.rm = TRUE), 
    sprintf("(%.1f%%)", 100 * mean(sam$WGD, na.rm = TRUE)), "\n")
cat("\nTop WGD frequencies:\n")
ct_wgd <- sam[, .(N = .N, WGD_Pct = round(100 * mean(WGD, na.rm = TRUE), 1)), 
              by = CANCER_TYPE][order(-WGD_Pct)]
print(head(ct_wgd, 10))
cat("\nTP53 mutation in WGD+:", 
    sprintf("%.1f%%", 100 * mean(sam$MUT_TP53[sam$WGD == 1])), "\n")
cat("TP53 mutation in WGD-:", 
    sprintf("%.1f%%", 100 * mean(sam$MUT_TP53[sam$WGD == 0])), "\n")

cat("\n=== Phase 2 Complete ===\n")
