# ============================================================
# Phase 1: Data Loading, QC, and Pan-Cancer WGD Frequency Atlas
# Project: Whole-Genome Doubling and Cancer Evolution
# Target: Cancer Cell Research Article
# ============================================================

.libPaths(c("/home/caiwj2001/R/library", .libPaths()))
suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(scales)
  library(ggrepel)
  library(gridExtra)
})

DATA_DIR <- "/home/caiwj2001/肿瘤突变特征与治疗反应/msk_impact_50k_2026"
OUT_DIR  <- "/home/caiwj2001/癌症演化"

# ---- 1. Load sample-level clinical data ----
cat("\n=== Loading sample data ===\n")
sam <- fread(file.path(DATA_DIR, "data_clinical_sample.txt"), skip = 4)
cat("Samples:", nrow(sam), "\n")
cat("Columns:", paste(names(sam), collapse=", "), "\n")

# ---- 2. Load patient-level data ----
pat <- fread(file.path(DATA_DIR, "data_clinical_patient.txt"), skip = 4)
cat("Patients:", nrow(pat), "\n")

# ---- 3. Merge patient data into samples ----
sam <- merge(sam, pat[, .(PATIENT_ID, SEX, OS_STATUS, OS_MONTHS, AGE_AT_DX, ANCESTRY_LABEL)], 
             by = "PATIENT_ID", all.x = TRUE)
cat("Merged samples:", nrow(sam), "\n")

# ---- 4. Data QC ----
cat("\n=== QC Report ===\n")
cat("Total samples:", nrow(sam), "\n")
cat("WGD=TRUE:", sum(sam$FACETS_WGD == "TRUE", na.rm=TRUE), 
    sprintf("(%.1f%%)", 100*mean(sam$FACETS_WGD == "TRUE", na.rm=TRUE)), "\n")
cat("WGD=FALSE:", sum(sam$FACETS_WGD == "FALSE", na.rm=TRUE), "\n")
cat("WGD=NA:", sum(is.na(sam$FACETS_WGD)), "\n\n")

# Sample types
cat("Sample types:\n")
print(table(sam$SAMPLE_TYPE, useNA="always"))

# FACETS QC
cat("\nFACETS_QC:\n")
print(table(sam$FACETS_QC, useNA="always"))

# Age QC
sam$AGE_AT_DX <- as.numeric(as.character(sam$AGE_AT_DX))
cat("\nAge summary:", summary(sam$AGE_AT_DX), "\n")

# Clean up OS
sam$OS_MONTHS <- as.numeric(sam$OS_MONTHS)
sam$OS_EVENT <- ifelse(sam$OS_STATUS == "DECEASED", 1, 0)
cat("\nOS events:", sum(sam$OS_EVENT, na.rm=TRUE), "/", sum(!is.na(sam$OS_EVENT)), "\n")

# ---- 5. Define WGD binary variable ----
sam$WGD <- ifelse(sam$FACETS_WGD == "TRUE", 1, 0)
sam$WGD <- ifelse(is.na(sam$FACETS_WGD), NA, sam$WGD)

# ---- 6. Pan-Cancer WGD Frequency ----
ct_wgd <- sam[, .(
  N_Total = .N,
  N_WGD = sum(WGD, na.rm = TRUE),
  WGD_Freq = round(100 * mean(WGD, na.rm = TRUE), 1)
), by = CANCER_TYPE]

ct_wgd <- ct_wgd[order(-N_Total)]
ct_wgd$SE <- round(sqrt(ct_wgd$WGD_Freq * (100 - ct_wgd$WGD_Freq) / ct_wgd$N_Total), 2)
ct_wgd$CI_Lower <- round(ct_wgd$WGD_Freq - 1.96 * ct_wgd$SE, 1)
ct_wgd$CI_Upper <- round(ct_wgd$WGD_Freq + 1.96 * ct_wgd$SE, 1)

cat("\n=== Pan-Cancer WGD Frequencies (by cancer type) ===\n")
print(ct_wgd, nrows = 40)

# Save
write.csv(ct_wgd, file.path(OUT_DIR, "tables", "TableS1_WGD_Frequency_by_CancerType.csv"), row.names = FALSE)

# ---- 7. WGD Frequency by Sample Type ----
cat("\n=== WGD by Sample Type ===\n")
st_wgd <- sam[, .(
  N = .N,
  WGD_N = sum(WGD, na.rm = TRUE),
  WGD_Pct = round(100 * mean(WGD, na.rm = TRUE), 1)
), by = SAMPLE_TYPE]
print(st_wgd)

# ---- 8. WGD Frequency by Genetic Ancestry ----
cat("\n=== WGD by Ancestry ===\n")
anc_wgd <- sam[, .(
  N = .N,
  WGD_Pct = round(100 * mean(WGD, na.rm = TRUE), 1)
), by = ANCESTRY_LABEL][N >= 100][order(-N)]
print(anc_wgd)

# ---- 9. Ploidy by WGD status ----
cat("\n=== Ploidy by WGD ===\n")
sam$FACETS_PLOIDY <- as.numeric(sam$FACETS_PLOIDY)
ploidy_wgd <- sam[!is.na(FACETS_PLOIDY) & !is.na(WGD), 
  .(N = .N, Mean_Ploidy = round(mean(FACETS_PLOIDY), 2),
    Median_Ploidy = median(FACETS_PLOIDY),
    SD_Ploidy = round(sd(FACETS_PLOIDY), 2)),
  by = WGD]
print(ploidy_wgd)

# ---- 10. TMB by WGD ----
sam$TMB_SCORE <- as.numeric(sam$TMB_SCORE)
cat("\n=== TMB by WGD (log10) ===\n")
tmb_wgd <- sam[!is.na(TMB_SCORE) & TMB_SCORE > 0 & !is.na(WGD), 
  .(N = .N, 
    Median_TMB = round(median(TMB_SCORE), 2),
    Mean_TMB = round(mean(TMB_SCORE), 2)),
  by = WGD]
print(tmb_wgd)

# Wilcoxon test for TMB ~ WGD
wtest <- wilcox.test(log10(TMB_SCORE) ~ WGD, data = sam[!is.na(WGD) & TMB_SCORE > 0])
cat("Wilcoxon TMB ~ WGD: p =", format(wtest$p.value, digits=3, scientific=TRUE), "\n")

# ---- 11. Save cleaned merged data for downstream analysis ----
saveRDS(sam, file.path(OUT_DIR, "data", "sample_merged.rds"))
cat("\nSaved merged sample data to", file.path(OUT_DIR, "data", "sample_merged.rds"), "\n")

cat("\n=== Phase 1 Complete ===\n")
