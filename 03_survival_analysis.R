# ============================================================
# Phase 3: Survival Analysis
# WGD Prognostic Impact: KM + Cox Regression
# Pan-Cancer + Cancer-Type-Specific
# ============================================================

.libPaths(c("/home/caiwj2001/R/library", .libPaths()))
suppressMessages({
  library(data.table)
  library(dplyr)
  library(survival)
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(grid)
})

OUT_DIR <- "/home/caiwj2001/癌症演化"

# Load expanded dataset
sam <- readRDS(file.path(OUT_DIR, "data", "sample_with_genomics.rds"))
cat("Loaded", nrow(sam), "samples\n")

# ---- 1. Build Survival Analysis Cohort ----
cat("\n=== Building survival cohort ===\n")

# Filter: Primary samples, complete OS data, known WGD
surv <- sam[SAMPLE_TYPE == "Primary" & !is.na(OS_MONTHS) & !is.na(OS_STATUS) & !is.na(WGD)]
cat("Primary samples with OS data:", nrow(surv), "\n")

# Clean up
surv$OS_MONTHS <- as.numeric(surv$OS_MONTHS)
surv$OS_EVENT <- ifelse(surv$OS_STATUS == "DECEASED", 1, 0)
surv$AGE_AT_DX <- as.numeric(as.character(surv$AGE_AT_DX))
surv$TMB_SCORE <- as.numeric(surv$TMB_SCORE)

# Remove NA or 0 survival times
surv <- surv[OS_MONTHS > 0]
cat("After removing OS_MONTHS <= 0:", nrow(surv), "\n")

# Remove NA age
surv <- surv[!is.na(AGE_AT_DX) & AGE_AT_DX > 0 & AGE_AT_DX < 120]
cat("After age filter:", nrow(surv), "\n")

# WGD factor
surv$WGD_factor <- factor(surv$WGD, levels = c(0, 1), labels = c("WGD-", "WGD+"))

cat("Final survival cohort:", nrow(surv), "\n")
cat("Events (deceased):", sum(surv$OS_EVENT), "\n")
cat("WGD+:", sum(surv$WGD == 1), "(", round(100*mean(surv$WGD), 1), "%)\n")

# ---- 2. Pan-Cancer Kaplan-Meier ----
cat("\n=== Pan-Cancer Kaplan-Meier ===\n")
fit_km <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ WGD_factor, data = surv)
print(fit_km)

# Log-rank test
lr_test <- survdiff(Surv(OS_MONTHS, OS_EVENT) ~ WGD_factor, data = surv)
cat("\nLog-rank test:\n")
print(lr_test)

# Median survival
cat("\nMedian OS:\n")
print(summary(fit_km)$table[, c("records", "events", "median")])

# ---- 3. Multivariate Cox Regression (pan-cancer) ----
cat("\n=== Multivariate Cox Regression ===\n")

# Log-transform TMB
surv$log10TMB <- log10(surv$TMB_SCORE + 0.01)

# Model 1: Unadjusted
cox1 <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ WGD, data = surv)
cat("\nModel 1 (unadjusted):\n")
print(summary(cox1)$coefficients)

# Model 2: Adjusted for age + sex
cox2 <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ WGD + AGE_AT_DX + SEX, data = surv)
cat("\nModel 2 (age + sex adjusted):\n")
print(summary(cox2)$coefficients)

# Model 3: Fully adjusted (age, sex, TMB, TP53)
cox3 <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ WGD + AGE_AT_DX + SEX + log10TMB + MUT_TP53, data = surv)
cat("\nModel 3 (fully adjusted):\n")
print(summary(cox3)$coefficients)

# ---- 4. Cancer-Type-Specific Cox Regression ----
cat("\n=== Cancer-Type-Specific Cox Regression ===\n")

# Cancer types with >= 100 primary samples and >= 20 events
ct_counts <- surv[, .(N = .N, Events = sum(OS_EVENT)), by = CANCER_TYPE]
ct_valid <- ct_counts[N >= 100 & Events >= 20, CANCER_TYPE]
cat("Cancer types with >=100 samples and >=20 events:", length(ct_valid), "\n")

ct_results <- data.frame()

for (ct in ct_valid) {
  sub <- surv[CANCER_TYPE == ct & !is.na(AGE_AT_DX)]
  wgd_plus <- sum(sub$WGD == 1)
  wgd_minus <- sum(sub$WGD == 0)
  
  if (wgd_plus >= 10 && wgd_minus >= 10) {
    # Unadjusted
    fit_u <- tryCatch({
      coxph(Surv(OS_MONTHS, OS_EVENT) ~ WGD, data = sub)
    }, error = function(e) NULL)
    
    if (!is.null(fit_u)) {
      coef_u <- summary(fit_u)$coefficients
      wgd_row <- which(rownames(coef_u) == "WGD")
      if (length(wgd_row) == 1) {
        ct_results <- rbind(ct_results, data.frame(
          Cancer_Type = ct,
          N = nrow(sub),
          Events = sum(sub$OS_EVENT),
          WGD_plus_N = wgd_plus,
          WGD_plus_Pct = round(100 * wgd_plus / nrow(sub), 1),
          HR_unadj = round(coef_u[wgd_row, "exp(coef)"], 3),
          HR_unadj_lower = round(exp(coef_u[wgd_row, "coef"] - 1.96 * coef_u[wgd_row, "se(coef)"]), 3),
          HR_unadj_upper = round(exp(coef_u[wgd_row, "coef"] + 1.96 * coef_u[wgd_row, "se(coef)"]), 3),
          P_unadj = coef_u[wgd_row, "Pr(>|z|)"],
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

ct_results <- ct_results[order(ct_results$HR_unadj), ]
ct_results$FDR <- p.adjust(ct_results$P_unadj, method = "BH")
cat("\nCancer-type-specific HRs:\n")
print(ct_results, nrows = 30)

# Save
write.csv(ct_results, file.path(OUT_DIR, "tables", "TableS3_CancerType_Cox_Regression.csv"), row.names = FALSE)

# ---- 5. KM by TP53/WGD combined ----
cat("\n=== KM by TP53/WGD combined ===\n")
surv$group_TP53_WGD <- paste0(
  ifelse(surv$MUT_TP53 == 1, "TP53mut", "TP53wt"),
  "_",
  ifelse(surv$WGD == 1, "WGD+", "WGD-")
)
surv$group_TP53_WGD <- factor(surv$group_TP53_WGD, 
  levels = c("TP53wt_WGD-", "TP53wt_WGD+", "TP53mut_WGD-", "TP53mut_WGD+"))

fit_k4 <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ group_TP53_WGD, data = surv)
cat("\nTP53/WGD 4-group KM:\n")
print(fit_k4, print.rmean = TRUE)

# ---- 6. Cox for TP53 × WGD interaction ----
cox_int <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ WGD * MUT_TP53 + AGE_AT_DX + SEX + log10TMB, data = surv)
cat("\nInteraction model (WGD × TP53):\n")
print(summary(cox_int)$coefficients)

# ---- 7. WGD combined with TMB strata ----
surv$TMB_group <- cut(surv$TMB_SCORE, 
  breaks = c(-Inf, 5, 10, 20, Inf),
  labels = c("TMB<5", "TMB 5-10", "TMB 10-20", "TMB>20"))

surv$WGD_TMB_group <- paste0(
  ifelse(surv$WGD == 1, "WGD+", "WGD-"),
  " / ",
  surv$TMB_group
)

fit_wgdtmb <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ WGD_TMB_group, data = surv)
cat("\nWGD/TMB stratification:\n")
print(summary(fit_wgdtmb)$table[, c("records", "events", "median")])

# ---- 8. Save survival cohort ----
saveRDS(surv, file.path(OUT_DIR, "data", "survival_cohort.rds"))
cat("\nSaved survival cohort:", nrow(surv), "samples\n")

# ---- 9. Key results for manuscript ----
cat("\n========================================\n")
cat("SURVIVAL KEY RESULTS\n")
cat("========================================\n")
cat("Pan-cancer WGD HR (unadjusted):", 
    round(summary(cox1)$coefficients["WGD", "exp(coef)"], 3),
    sprintf("(%.3f-%.3f, P=%.2e)", 
            summary(cox1)$conf.int["WGD", "lower .95"],
            summary(cox1)$conf.int["WGD", "upper .95"],
            summary(cox1)$coefficients["WGD", "Pr(>|z|)"]), "\n")
cat("Pan-cancer WGD HR (fully adjusted):", 
    round(summary(cox3)$coefficients["WGD", "exp(coef)"], 3),
    sprintf("(%.3f-%.3f, P=%.2e)", 
            summary(cox3)$conf.int["WGD", "lower .95"],
            summary(cox3)$conf.int["WGD", "upper .95"],
            summary(cox3)$coefficients["WGD", "Pr(>|z|)"]), "\n")

# Significant cancer types
# Significant cancer types (using FDR-adjusted)
sig_ct <- ct_results[p.adjust(ct_results$P_unadj, method = "BH") < 0.05, ]
cat("\nCancer types with significant WGD prognostic effect (FDR<0.05):\n")
if (nrow(sig_ct) > 0) {
  sig_ct$FDR_val <- p.adjust(sig_ct$P_unadj, method = "BH")
  print(sig_ct)
} else {
  cat("None at FDR<0.05\n")
}
cat("\nAll CT results with nominal P<0.05:\n")
nom_sig <- ct_results[ct_results$P_unadj < 0.05, ]
print(nom_sig)

cat("\n=== Phase 3 Complete ===\n")
