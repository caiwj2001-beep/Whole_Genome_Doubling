# ============================================================
# Phase 4: Comprehensive Figure Generation
# 8 Publication-Quality Figures (TIFF 600dpi + PNG 300dpi + PDF)
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
  library(grid)
})

OUT_DIR <- "/home/caiwj2001/癌症演化"
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Load data
sam <- readRDS(file.path(OUT_DIR, "data", "sample_with_genomics.rds"))
surv <- readRDS(file.path(OUT_DIR, "data", "survival_cohort.rds"))
cat("Data loaded\n")

# Common theme
theme_cell <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

# Color scheme
col_wgd_plus  <- "#E74C3C"
col_wgd_minus <- "#3498DB"
col_neutral   <- "#95A5A6"

# ============================================================
# Figure 1: Pan-Cancer WGD Frequency Atlas
# ============================================================
cat("\n=== Figure 1: Pan-Cancer WGD Frequency ===\n")

fig1_data <- sam[, .(
  N = .N,
  WGD_Pct = 100 * mean(WGD, na.rm = TRUE),
  WGD_N = sum(WGD, na.rm = TRUE)
), by = CANCER_TYPE][N >= 50][order(-WGD_Pct)]

fig1_data$CANCER_TYPE <- factor(fig1_data$CANCER_TYPE, 
  levels = fig1_data$CANCER_TYPE[order(fig1_data$WGD_Pct)])

fig1_data$WGD_SE <- sqrt(fig1_data$WGD_Pct * (100 - fig1_data$WGD_Pct) / fig1_data$N)

p1 <- ggplot(fig1_data, aes(x = WGD_Pct, y = CANCER_TYPE)) +
  geom_bar(stat = "identity", fill = col_wgd_plus, alpha = 0.85, width = 0.7) +
  geom_errorbarh(aes(xmin = pmax(0, WGD_Pct - 1.96 * WGD_SE), 
                     xmax = pmin(100, WGD_Pct + 1.96 * WGD_SE)),
                 height = 0.2, linewidth = 0.3) +
  geom_text(aes(label = paste0(sprintf("%.1f", WGD_Pct), "%  (", N, ")")),
            hjust = -0.1, size = 2.5, color = "grey30") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(fig1_data$WGD_Pct) * 1.35)) +
  labs(x = "WGD Frequency (%)", y = "", 
       title = "Pan-Cancer Whole-Genome Doubling Frequency") +
  theme_cell +
  theme(axis.text.y = element_text(size = 7.5))

ggsave(file.path(FIG_DIR, "Figure1_WGD_Frequency.pdf"), p1, width = 7.5, height = 8, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure1_WGD_Frequency.png"), p1, width = 7.5, height = 8, dpi = 300)
tiff(file.path(FIG_DIR, "Figure1_WGD_Frequency.tiff"), width = 7.5, height = 8, 
     units = "in", res = 600, compression = "lzw")
print(p1)
dev.off()
cat("Figure 1 saved\n")

# ============================================================
# Figure 2: Driver Gene Associations with WGD (Volcano Plot)
# ============================================================
cat("=== Figure 2: Driver Gene Associations ===\n")

assoc <- read.csv(file.path(OUT_DIR, "tables", "TableS2_WGD_DriverGene_Associations.csv"))
assoc$log2OR <- log2(assoc$OR)
assoc$negLog10FDR <- -log10(assoc$FDR)
assoc$negLog10FDR[is.infinite(assoc$negLog10FDR)] <- max(assoc$negLog10FDR[is.finite(assoc$negLog10FDR)]) + 2
assoc$direction <- ifelse(assoc$OR > 1, "Enriched in WGD+", "Depleted in WGD+")
assoc$label <- ifelse(assoc$FDR < 0.05, assoc$Gene, "")

p2 <- ggplot(assoc, aes(x = log2OR, y = negLog10FDR, color = direction, label = label)) +
  geom_point(aes(size = abs(log2OR)), alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 25, box.padding = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("Enriched in WGD+" = col_wgd_plus, 
                                 "Depleted in WGD+" = col_wgd_minus)) +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(x = "log2(Odds Ratio)", y = "-log10(FDR)", 
       title = "Driver Gene Alterations Associated with WGD",
       color = "") +
  theme_cell

ggsave(file.path(FIG_DIR, "Figure2_Driver_Associations.pdf"), p2, width = 7, height = 6, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure2_Driver_Associations.png"), p2, width = 7, height = 6, dpi = 300)
tiff(file.path(FIG_DIR, "Figure2_Driver_Associations.tiff"), width = 7, height = 6, 
     units = "in", res = 600, compression = "lzw")
print(p2)
dev.off()
cat("Figure 2 saved\n")

# ============================================================
# Figure 3: Ploidy & TMB by WGD Status
# ============================================================
cat("=== Figure 3: Ploidy and TMB ===\n")

sam_plot <- sam[!is.na(WGD) & !is.na(FACETS_PLOIDY) & FACETS_PLOIDY > 0 & FACETS_PLOIDY < 10]
sam_plot$WGD_label <- ifelse(sam_plot$WGD == 1, "WGD+", "WGD-")

p3a <- ggplot(sam_plot, aes(x = WGD_label, y = FACETS_PLOIDY, fill = WGD_label)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c("WGD+" = col_wgd_plus, "WGD-" = col_wgd_minus)) +
  labs(x = "", y = "Ploidy (FACETS)", title = "Ploidy by WGD Status") +
  theme_cell + theme(legend.position = "none")

sam_tmb <- sam[!is.na(WGD) & TMB_SCORE > 0 & TMB_SCORE < 100]
sam_tmb$WGD_label <- ifelse(sam_tmb$WGD == 1, "WGD+", "WGD-")
sam_tmb$log10TMB <- log10(sam_tmb$TMB_SCORE)

p3b <- ggplot(sam_tmb, aes(x = WGD_label, y = log10TMB, fill = WGD_label)) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c("WGD+" = col_wgd_plus, "WGD-" = col_wgd_minus)) +
  labs(x = "", y = "log10(TMB)", title = "TMB by WGD Status") +
  theme_cell + theme(legend.position = "none")

p3 <- grid.arrange(p3a, p3b, ncol = 2)
ggsave(file.path(FIG_DIR, "Figure3_Ploidy_TMB.pdf"), p3, width = 8, height = 4, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure3_Ploidy_TMB.png"), p3, width = 8, height = 4, dpi = 300)
tiff(file.path(FIG_DIR, "Figure3_Ploidy_TMB.tiff"), width = 8, height = 4, 
     units = "in", res = 600, compression = "lzw")
grid.draw(p3)
dev.off()
cat("Figure 3 saved\n")

# ============================================================
# Figure 4: TP53/RB1 Combined WGD Rates
# ============================================================
cat("=== Figure 4: TP53/RB1 Interaction ===\n")

tp53rb1_data <- sam[, .(
  WGD_Rate = 100 * mean(WGD, na.rm = TRUE),
  N = .N
), by = .(MUT_TP53, MUT_RB1)]

tp53rb1_data$Group <- paste0(
  "TP53", ifelse(tp53rb1_data$MUT_TP53, "mut", "wt"),
  "\nRB1", ifelse(tp53rb1_data$MUT_RB1, "mut", "wt")
)

tp53rb1_data <- tp53rb1_data[order(-tp53rb1_data$WGD_Rate)]

p4 <- ggplot(tp53rb1_data, aes(x = reorder(Group, -WGD_Rate), y = WGD_Rate)) +
  geom_bar(stat = "identity", aes(fill = WGD_Rate), width = 0.6) +
  geom_text(aes(label = paste0(sprintf("%.1f", WGD_Rate), "%\n(N=", N, ")")), 
            vjust = -0.2, size = 3.5) +
  scale_fill_gradient(low = col_wgd_minus, high = col_wgd_plus, guide = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55)) +
  labs(x = "", y = "WGD Frequency (%)", 
       title = "WGD Frequency by TP53 and RB1 Mutation Status") +
  theme_cell

ggsave(file.path(FIG_DIR, "Figure4_TP53_RB1_WGD.pdf"), p4, width = 5, height = 4.5, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure4_TP53_RB1_WGD.png"), p4, width = 5, height = 4.5, dpi = 300)
tiff(file.path(FIG_DIR, "Figure4_TP53_RB1_WGD.tiff"), width = 5, height = 4.5, 
     units = "in", res = 600, compression = "lzw")
print(p4)
dev.off()
cat("Figure 4 saved\n")

# ============================================================
# Figure 5: Kaplan-Meier Survival Curves
# ============================================================
cat("=== Figure 5: Kaplan-Meier Curves ===\n")

# 5a: Pan-cancer KM
fit_km <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ WGD_factor, data = surv)
km_data <- data.frame(
  time = fit_km$time,
  surv = fit_km$surv,
  upper = fit_km$upper,
  lower = fit_km$lower,
  strata = rep(names(fit_km$strata), fit_km$strata)
)
km_data$strata <- gsub("WGD_factor=", "", km_data$strata)

p5a <- ggplot(km_data, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linewidth = 0) +
  scale_color_manual(values = c("WGD-" = col_wgd_minus, "WGD+" = col_wgd_plus)) +
  scale_fill_manual(values = c("WGD-" = col_wgd_minus, "WGD+" = col_wgd_plus)) +
  scale_x_continuous(breaks = seq(0, 120, 24), limits = c(0, 120)) +
  labs(x = "Overall Survival (months)", y = "Survival Probability",
       title = "Pan-Cancer Overall Survival by WGD Status",
       subtitle = "29,777 primary tumors; Log-rank P < 2e-16",
       color = "", fill = "") +
  theme_cell +
  annotate("text", x = 80, y = 0.9, 
           label = paste0("WGD-: median ", round(summary(fit_km)$table[1, "median"], 1), " mo"),
           color = col_wgd_minus, size = 3, hjust = 0) +
  annotate("text", x = 80, y = 0.82, 
           label = paste0("WGD+: median ", round(summary(fit_km)$table[2, "median"], 1), " mo"),
           color = col_wgd_plus, size = 3, hjust = 0)

# 5b: TP53/WGD 4-group KM
fit_k4 <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ group_TP53_WGD, data = surv)
km4_data <- data.frame(
  time = fit_k4$time,
  surv = fit_k4$surv,
  strata = rep(names(fit_k4$strata), fit_k4$strata)
)
km4_data$strata <- gsub("group_TP53_WGD=", "", km4_data$strata)

grp_colors <- c("TP53wt_WGD-" = "#2ECC71", "TP53wt_WGD+" = "#E74C3C",
                "TP53mut_WGD-" = "#3498DB", "TP53mut_WGD+" = "#8E44AD")

p5b <- ggplot(km4_data, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 0.7) +
  scale_color_manual(values = grp_colors) +
  scale_x_continuous(breaks = seq(0, 120, 24), limits = c(0, 120)) +
  labs(x = "Overall Survival (months)", y = "Survival Probability",
       title = "Survival Stratified by TP53 and WGD Status",
       color = "") +
  theme_cell

p5 <- grid.arrange(p5a, p5b, ncol = 1, heights = c(1, 1))
ggsave(file.path(FIG_DIR, "Figure5_Kaplan_Meier.pdf"), p5, width = 7, height = 8, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure5_Kaplan_Meier.png"), p5, width = 7, height = 8, dpi = 300)
tiff(file.path(FIG_DIR, "Figure5_Kaplan_Meier.tiff"), width = 7, height = 8, 
     units = "in", res = 600, compression = "lzw")
grid.draw(p5)
dev.off()
cat("Figure 5 saved\n")

# ============================================================
# Figure 6: Cancer-Type-Specific Forest Plot
# ============================================================
cat("=== Figure 6: Forest Plot ===\n")

ct_results <- read.csv(file.path(OUT_DIR, "tables", "TableS3_CancerType_Cox_Regression.csv"))
ct_results <- ct_results[order(-ct_results$HR_unadj), ]
ct_results$Cancer_Type <- factor(ct_results$Cancer_Type, 
  levels = ct_results$Cancer_Type[order(ct_results$HR_unadj)])

# Only show cancer types with enough precision
ct_plot <- ct_results[ct_results$N >= 100, ]
ct_plot$Significant <- ifelse(ct_plot$P_unadj < 0.05, "P < 0.05", "n.s.")

p6 <- ggplot(ct_plot, aes(x = HR_unadj, y = Cancer_Type, color = Significant)) +
  geom_point(aes(size = N / 50), alpha = 0.8) +
  geom_errorbarh(aes(xmin = HR_unadj_lower, xmax = HR_unadj_upper), 
                 height = 0.2, linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("P < 0.05" = col_wgd_plus, "n.s." = col_neutral)) +
  scale_size_continuous(range = c(1.5, 6), guide = "none") +
  labs(x = "Hazard Ratio (WGD+ vs WGD-)", y = "",
       title = "Cancer-Type-Specific Prognostic Impact of WGD",
       subtitle = "Unadjusted Cox regression; point size ∝ sample size",
       color = "") +
  theme_cell +
  annotate("text", x = max(ct_plot$HR_unadj_upper) * 0.85, 
           y = 1, label = "← WGD+ favorable   |   WGD+ unfavorable →",
           size = 3, color = "grey50")

ggsave(file.path(FIG_DIR, "Figure6_Forest_Plot.pdf"), p6, width = 8, height = 7, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure6_Forest_Plot.png"), p6, width = 8, height = 7, dpi = 300)
tiff(file.path(FIG_DIR, "Figure6_Forest_Plot.tiff"), width = 8, height = 7, 
     units = "in", res = 600, compression = "lzw")
print(p6)
dev.off()
cat("Figure 6 saved\n")

# ============================================================
# Figure 7: WGD × TMB Combined Stratification
# ============================================================
cat("=== Figure 7: WGD × TMB Stratification ===\n")

surv$TMB_strata <- ifelse(surv$TMB_SCORE >= 10, "TMB-High (≥10)", "TMB-Low (<10)")
surv$WGD_TMB_combo <- paste0(
  ifelse(surv$WGD == 1, "WGD+", "WGD-"), " / ",
  surv$TMB_strata
)

fit_tmb <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ WGD_TMB_combo, data = surv)
tmb_km_data <- data.frame(
  time = fit_tmb$time,
  surv = fit_tmb$surv,
  strata = rep(names(fit_tmb$strata), fit_tmb$strata)
)
tmb_km_data$strata <- gsub("WGD_TMB_combo=", "", tmb_km_data$strata)

combo_colors <- c("WGD- / TMB-Low (<10)" = "#3498DB", 
                  "WGD+ / TMB-Low (<10)" = "#E74C3C",
                  "WGD- / TMB-High (≥10)" = "#1ABC9C",
                  "WGD+ / TMB-High (≥10)" = "#F39C12")

p7 <- ggplot(tmb_km_data, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 0.7) +
  scale_color_manual(values = combo_colors) +
  scale_x_continuous(breaks = seq(0, 120, 24), limits = c(0, 120)) +
  labs(x = "Overall Survival (months)", y = "Survival Probability",
       title = "Combined WGD and TMB Stratification",
       color = "") +
  theme_cell

ggsave(file.path(FIG_DIR, "Figure7_WGD_TMB_Stratification.pdf"), p7, width = 7, height = 5, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure7_WGD_TMB_Stratification.png"), p7, width = 7, height = 5, dpi = 300)
tiff(file.path(FIG_DIR, "Figure7_WGD_TMB_Stratification.tiff"), width = 7, height = 5, 
     units = "in", res = 600, compression = "lzw")
print(p7)
dev.off()
cat("Figure 7 saved\n")

# ============================================================
# Figure 8: WGD by Metastatic vs Primary
# ============================================================
cat("=== Figure 8: Metastatic Enrichment ===\n")

met_data <- sam[SAMPLE_TYPE %in% c("Primary", "Metastasis") & !is.na(WGD)]
met_summary <- met_data[, .(
  WGD_Pct = 100 * mean(WGD, na.rm = TRUE),
  N = .N,
  WGD_N = sum(WGD)
), by = .(SAMPLE_TYPE, CANCER_TYPE)]

# Top cancer types by sample count
top_ct <- sam[, .N, by = CANCER_TYPE][order(-N)][1:15, CANCER_TYPE]
met_plot <- met_summary[CANCER_TYPE %in% top_ct & N >= 20]

met_plot$CANCER_TYPE <- factor(met_plot$CANCER_TYPE, 
  levels = top_ct[top_ct %in% unique(met_plot$CANCER_TYPE)])

p8 <- ggplot(met_plot, aes(x = CANCER_TYPE, y = WGD_Pct, fill = SAMPLE_TYPE)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c("Primary" = col_wgd_minus, "Metastasis" = col_wgd_plus)) +
  labs(x = "", y = "WGD Frequency (%)", fill = "",
       title = "WGD Enrichment in Metastatic vs Primary Tumors") +
  theme_cell +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(FIG_DIR, "Figure8_Metastatic_Enrichment.pdf"), p8, width = 8, height = 5, device = "pdf")
ggsave(file.path(FIG_DIR, "Figure8_Metastatic_Enrichment.png"), p8, width = 8, height = 5, dpi = 300)
tiff(file.path(FIG_DIR, "Figure8_Metastatic_Enrichment.tiff"), width = 8, height = 5, 
     units = "in", res = 600, compression = "lzw")
print(p8)
dev.off()
cat("Figure 8 saved\n")

# ---- Verify all figures ----
cat("\n=== Figure verification ===\n")
fig_files <- list.files(FIG_DIR, pattern = "\\.(tiff|png|pdf)$")
for (f in sort(fig_files)) {
  info <- file.info(file.path(FIG_DIR, f))
  cat(sprintf("%-40s %8.0f bytes\n", f, info$size))
}

cat("\n=== Phase 4 Complete ===\n")
