#!/usr/bin/env Rscript

# 02n_model_validation_sensitivity.R
# Additional model validation and sensitivity analyses:
# - 3-year landmark calibration intercept/slope for transported risks
# - 3-year external DCA and fixed-threshold net benefit
# - Bootstrap feature-selection stability in the discovery cohort
# - TCGA proliferation baseline/sensitivity analyses
# - BH-FDR columns for multi-omics and immune summary tables

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(glmnet)
  library(ggplot2)
  library(hgu133a2.db)
  library(AnnotationDbi)
})

proc_dir <- "data/processed"
ref_dir <- "data/references"
res_dir <- "results"
plot_dir <- "plots/supplementary"
bench_dir <- file.path(res_dir, "benchmarks")
effect_dir <- file.path(res_dir, "effect_sizes")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bench_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(effect_dir, showWarnings = FALSE, recursive = TRUE)

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_risk_csv <- function(path, dataset) {
  df <- safe_read_csv(path)
  if ("time" %in% names(df) && !"time_months" %in% names(df)) df$time_months <- df$time
  if (!all(c("risk_score", "time_months", "status") %in% names(df))) {
    stop("Risk file lacks risk_score/time_months/status columns: ", path)
  }
  df %>%
    transmute(
      dataset = dataset,
      id = if ("sample" %in% names(df)) sample else if ("patient_id" %in% names(df)) patient_id else row_number(),
      time_months = as.numeric(time_months),
      status = as.integer(status),
      risk_score = as.numeric(risk_score)
    ) %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))
}

landmark_df <- function(df, t_months = 36) {
  df %>%
    mutate(
      eligible_3y = (time_months > t_months) | (time_months <= t_months & status == 1),
      event_3y = as.integer(time_months <= t_months & status == 1)
    ) %>%
    filter(eligible_3y)
}

ipcw_brier <- function(time_months, status, pred_risk, t_months) {
  df <- tibble(
    time_months = as.numeric(time_months),
    status = as.integer(status),
    pred_risk = as.numeric(pred_risk)
  ) %>%
    filter(!is.na(time_months), !is.na(status), !is.na(pred_risk))

  gfit <- survfit(Surv(df$time_months, 1 - df$status) ~ 1)
  gstep <- approxfun(gfit$time, gfit$surv, method = "constant", f = 0, rule = 2)
  y <- as.integer(df$time_months <= t_months & df$status == 1)
  w <- rep(0, nrow(df))
  is_event_before <- df$time_months <= t_months & df$status == 1
  is_at_risk_t <- df$time_months > t_months
  w[is_event_before] <- 1 / pmax(as.numeric(gstep(df$time_months[is_event_before])), 1e-12)
  w[is_at_risk_t] <- 1 / pmax(as.numeric(gstep(t_months)), 1e-12)
  mean(w * (y - df$pred_risk) ^ 2, na.rm = TRUE)
}

calc_nb <- function(pred_risk, outcome, threshold) {
  pred_pos <- pred_risk >= threshold
  n <- length(outcome)
  tp <- sum(pred_pos & outcome == 1, na.rm = TRUE)
  fp <- sum(pred_pos & outcome == 0, na.rm = TRUE)
  (tp / n) - (fp / n) * threshold / (1 - threshold)
}

bootstrap_nb <- function(pred_risk, outcome, threshold, n_boot = 300, seed = 20260628) {
  set.seed(seed)
  n <- length(outcome)
  vals <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    vals[i] <- calc_nb(pred_risk[idx], outcome[idx], threshold)
  }
  qs <- quantile(vals, c(0.025, 0.975), na.rm = TRUE)
  c(mean = mean(vals, na.rm = TRUE), ci_low = unname(qs[[1]]), ci_high = unname(qs[[2]]))
}

zscore_by_gene <- function(expr_mat) {
  z <- t(scale(t(expr_mat)))
  z[is.na(z)] <- 0
  z
}

build_gse14520_gene_expr <- function() {
  expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
  gene_symbols <- mapIds(
    hgu133a2.db,
    keys = rownames(expr_14),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  keep <- !is.na(gene_symbols) & gene_symbols != ""
  expr_raw <- expr_14[keep, , drop = FALSE]
  gene_symbols <- gene_symbols[keep]
  expr_sum <- rowsum(expr_raw, group = gene_symbols, reorder = FALSE)
  expr_cnt <- rowsum(matrix(1, nrow(expr_raw), ncol(expr_raw)), group = gene_symbols, reorder = FALSE)
  expr_sum / expr_cnt
}

message("[Model validation] 3-year transported-risk calibration and decision-curve analysis...")

train_df <- read_risk_csv(file.path(res_dir, "risk_score_data.csv"), "GSE14520")
cox_offset <- coxph(Surv(time_months, status) ~ offset(risk_score), data = train_df, x = TRUE, y = TRUE)
bh <- basehaz(cox_offset, centered = FALSE)
H0_step <- approxfun(bh$time, bh$hazard, method = "constant", f = 0, rule = 2)
predict_risk <- function(lp, t_months) 1 - exp(-as.numeric(H0_step(t_months)) * exp(lp))

cohort_files <- c(
  "GSE14520" = file.path(res_dir, "risk_score_data.csv"),
  "TCGA-LIHC" = file.path(res_dir, "TCGA_LIHC_risk_score.csv"),
  "GSE76427" = file.path(res_dir, "GSE76427_risk_score.csv"),
  "GSE10143-HCC" = file.path(res_dir, "GSE10143_HCC_risk_score.csv"),
  "GSE27150" = file.path(res_dir, "GSE27150_risk_score.csv"),
  "ICGC-LIRI-JP (HCCDB18)" = file.path(res_dir, "HCCDB18_LIRIJP_risk_score.csv")
)

t_landmark <- 36
risk_frames <- lapply(names(cohort_files), function(nm) {
  df <- read_risk_csv(cohort_files[[nm]], nm)
  df$pred_risk_3y <- pmin(pmax(predict_risk(df$risk_score, t_landmark), 1e-6), 1 - 1e-6)
  df
})
names(risk_frames) <- names(cohort_files)

cal_rows <- lapply(names(risk_frames), function(nm) {
  ldf <- landmark_df(risk_frames[[nm]], t_landmark)
  ldf$logit_pred <- qlogis(ldf$pred_risk_3y)
  out <- tibble(
    Dataset = nm,
    N_total = nrow(risk_frames[[nm]]),
    N_landmark = nrow(ldf),
    Events_3y = sum(ldf$event_3y == 1),
    Predicted_risk_3y_mean = mean(ldf$pred_risk_3y),
    Observed_event_3y = mean(ldf$event_3y),
    Calibration_intercept = NA_real_,
    Calibration_intercept_SE = NA_real_,
    Calibration_slope = NA_real_,
    Calibration_slope_SE = NA_real_,
    Calibration_slope_P = NA_real_,
    Brier_3y_landmark = mean((ldf$event_3y - ldf$pred_risk_3y) ^ 2),
    Brier_3y_IPCW = ipcw_brier(risk_frames[[nm]]$time_months, risk_frames[[nm]]$status, risk_frames[[nm]]$pred_risk_3y, t_landmark)
  )
  if (nrow(ldf) >= 20 && length(unique(ldf$event_3y)) == 2) {
    intercept_fit <- suppressWarnings(glm(event_3y ~ offset(logit_pred), data = ldf, family = binomial()))
    slope_fit <- suppressWarnings(glm(event_3y ~ logit_pred, data = ldf, family = binomial()))
    out$Calibration_intercept <- coef(intercept_fit)[[1]]
    out$Calibration_intercept_SE <- summary(intercept_fit)$coefficients[1, 2]
    out$Calibration_slope <- coef(slope_fit)[["logit_pred"]]
    out$Calibration_slope_SE <- summary(slope_fit)$coefficients["logit_pred", 2]
    out$Calibration_slope_P <- summary(slope_fit)$coefficients["logit_pred", 4]
  }
  out
})
calibration_df <- bind_rows(cal_rows)
write.csv(calibration_df, file.path(res_dir, "external_calibration_slope_intercept_3y.csv"), row.names = FALSE)

fixed_threshold <- median(landmark_df(risk_frames[["GSE14520"]], t_landmark)$pred_risk_3y, na.rm = TRUE)
threshold_grid <- seq(0.01, 0.60, by = 0.01)

dca_rows <- lapply(names(risk_frames), function(nm) {
  ldf <- landmark_df(risk_frames[[nm]], t_landmark)
  bind_rows(lapply(threshold_grid, function(th) {
    event_rate <- mean(ldf$event_3y)
    tibble(
      Dataset = nm,
      Threshold = th,
      NetBenefit_Model = calc_nb(ldf$pred_risk_3y, ldf$event_3y, th),
      NetBenefit_All = event_rate - (1 - event_rate) * th / (1 - th),
      NetBenefit_None = 0
    )
  }))
})
dca_df <- bind_rows(dca_rows)
write.csv(dca_df, file.path(res_dir, "external_dca_3y.csv"), row.names = FALSE)

threshold_rows <- lapply(names(risk_frames), function(nm) {
  ldf <- landmark_df(risk_frames[[nm]], t_landmark)
  event_rate <- mean(ldf$event_3y)
  nb <- bootstrap_nb(ldf$pred_risk_3y, ldf$event_3y, fixed_threshold)
  tibble(
    Dataset = nm,
    N_landmark = nrow(ldf),
    Events_3y = sum(ldf$event_3y == 1),
    Fixed_threshold_source = "Discovery median transported 3-year risk",
    Fixed_threshold = fixed_threshold,
    NetBenefit_Model = calc_nb(ldf$pred_risk_3y, ldf$event_3y, fixed_threshold),
    NetBenefit_Model_boot_mean = nb[["mean"]],
    NetBenefit_Model_CI95_L = nb[["ci_low"]],
    NetBenefit_Model_CI95_U = nb[["ci_high"]],
    NetBenefit_All = event_rate - (1 - event_rate) * fixed_threshold / (1 - fixed_threshold),
    NetBenefit_None = 0
  )
})
threshold_df <- bind_rows(threshold_rows)
write.csv(threshold_df, file.path(res_dir, "fixed_threshold_net_benefit_3y.csv"), row.names = FALSE)

dca_long <- dca_df %>%
  pivot_longer(starts_with("NetBenefit_"), names_to = "Strategy", values_to = "NetBenefit") %>%
  mutate(Strategy = recode(Strategy, NetBenefit_Model = "Model-guided", NetBenefit_All = "Intervene all", NetBenefit_None = "Intervene none"))

p_dca <- ggplot(dca_long, aes(Threshold, NetBenefit, color = Strategy, linewidth = Strategy)) +
  geom_line() +
  facet_wrap(~ Dataset, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("Model-guided" = "#2F6B7C", "Intervene all" = "#9B6A44", "Intervene none" = "#606060")) +
  scale_linewidth_manual(values = c("Model-guided" = 0.9, "Intervene all" = 0.55, "Intervene none" = 0.55)) +
  coord_cartesian(xlim = c(0, 0.60)) +
  labs(x = "Threshold probability", y = "Net benefit", color = NULL, linewidth = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(plot_dir, "Supp_Figure_S11_External_DCA_3y.pdf"), p_dca, width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(plot_dir, "Supp_Figure_S11_External_DCA_3y.png"), p_dca, width = 11, height = 7, dpi = 300, bg = "white")

message("[Model validation] Bootstrap feature-selection stability...")

expr_gene <- build_gse14520_gene_expr()
clinical_14 <- readRDS(file.path(proc_dir, "GSE14520_tumor_clinical.rds"))
matched_gsm <- intersect(colnames(expr_gene), clinical_14$Affy_GSM)
expr_matched <- expr_gene[, matched_gsm, drop = FALSE]
clinical_matched <- clinical_14[match(matched_gsm, clinical_14$Affy_GSM), ]
expr_scaled <- zscore_by_gene(expr_matched)

surv_data <- data.frame(
  sample = matched_gsm,
  time = as.numeric(clinical_matched$Survival.months),
  status = as.numeric(clinical_matched$Survival.status),
  stringsAsFactors = FALSE
) %>% filter(!is.na(time), !is.na(status))

deg_all <- safe_read_csv(file.path(res_dir, "deg_GSE14520_all.csv"))
ferro_path <- if (file.exists(file.path(res_dir, "ferroptosis_genes_hcc_context.csv"))) {
  file.path(res_dir, "ferroptosis_genes_hcc_context.csv")
} else {
  file.path(ref_dir, "ferroptosis_genes_expanded.csv")
}
candidate_genes <- intersect(
  unique(toupper(deg_all$Gene[deg_all$adj.P.Val < 0.05 & abs(deg_all$logFC) > 1])),
  unique(toupper(safe_read_csv(ferro_path)$Gene))
)
candidate_genes <- intersect(candidate_genes, rownames(expr_scaled))

set.seed(20260628)
n_boot_select <- suppressWarnings(as.integer(Sys.getenv("SELECTION_BOOTSTRAP_N", "200")))
if (is.na(n_boot_select) || n_boot_select < 50) n_boot_select <- 200

stability_path <- file.path(res_dir, "lasso_selection_stability.csv")
if (file.exists(stability_path)) {
  message("[Model validation] Reusing existing feature-selection stability table: ", stability_path)
  stability_df <- safe_read_csv(stability_path)
} else {
  selected_list <- vector("list", n_boot_select)
  for (i in seq_len(n_boot_select)) {
    if (i %% 20 == 0) message("[Model validation] Selection bootstrap ", i, "/", n_boot_select)
    idx <- sample(seq_len(nrow(surv_data)), nrow(surv_data), replace = TRUE)
    boot_surv <- surv_data[idx, , drop = FALSE]
    boot_samples <- boot_surv$sample
    cox_p <- rep(NA_real_, length(candidate_genes))
    names(cox_p) <- candidate_genes
    for (g in candidate_genes) {
      gx <- as.numeric(expr_scaled[g, boot_samples])
      if (sd(gx, na.rm = TRUE) <= 0.1) next
      fit <- try(coxph(Surv(time, status) ~ gx, data = boot_surv), silent = TRUE)
      if (!inherits(fit, "try-error")) cox_p[g] <- summary(fit)$coefficients[1, 5]
    }
    sig_genes <- names(cox_p)[is.finite(cox_p) & cox_p < 0.1]
    if (length(sig_genes) < 2) sig_genes <- names(sort(cox_p, na.last = NA))[seq_len(min(5, sum(is.finite(cox_p))))]
    if (length(sig_genes) >= 2) {
      x <- t(expr_scaled[sig_genes, boot_samples, drop = FALSE])
      y <- Surv(boot_surv$time, boot_surv$status)
      fit <- try(cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 5), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        cf <- coef(fit, s = fit$lambda.min)
        selected_list[[i]] <- rownames(cf)[which(as.numeric(cf) != 0)]
      }
    }
  }

  coef_df <- safe_read_csv(file.path(res_dir, "prognostic_model_coef.csv"))
  stability_df <- tibble(Gene = sort(unique(c(candidate_genes, coef_df$Gene)))) %>%
    mutate(
      Selected_N = vapply(Gene, function(g) sum(vapply(selected_list, function(x) g %in% x, logical(1))), integer(1)),
      Bootstrap_N = n_boot_select,
      Selection_proportion = Selected_N / Bootstrap_N,
      In_final_model = Gene %in% coef_df$Gene
    ) %>%
    arrange(desc(In_final_model), desc(Selection_proportion), Gene)
  write.csv(stability_df, stability_path, row.names = FALSE)
}

p_stability <- stability_df %>%
  filter(In_final_model | Selection_proportion >= 0.10) %>%
  mutate(Gene = reorder(Gene, Selection_proportion)) %>%
  ggplot(aes(Gene, Selection_proportion, fill = In_final_model)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "#B7C2C8", "TRUE" = "#2F6B7C")) +
  labs(x = NULL, y = "Bootstrap selection proportion", fill = "Final model") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(plot_dir, "Supp_Figure_S12_LASSO_Selection_Stability.pdf"), p_stability, width = 7.5, height = 6, device = cairo_pdf)
ggsave(file.path(plot_dir, "Supp_Figure_S12_LASSO_Selection_Stability.png"), p_stability, width = 7.5, height = 6, dpi = 300, bg = "white")

message("[Model validation] TCGA proliferation baseline and adjustment...")

tcga_risk <- read_risk_csv(file.path(res_dir, "TCGA_LIHC_risk_score.csv"), "TCGA-LIHC")
expr_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_expr.rds"))
patient_id <- substr(colnames(expr_tcga), 1, 12)
expr_tcga_t <- t(expr_tcga)
expr_tcga_sum <- rowsum(expr_tcga_t, patient_id)
expr_tcga_avg <- t(expr_tcga_sum / as.vector(table(patient_id)))

extract_gene_z <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(expr_mat)))
  z <- zscore_by_gene(expr_mat[genes, , drop = FALSE])
  colMeans(z, na.rm = TRUE)
}

cell_cycle_genes <- c("MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "CDC20", "BIRC5", "AURKA", "AURKB", "CCNA2", "MCM2", "MCM6")
tcga_prolif <- tibble(
  id = colnames(expr_tcga_avg),
  MKI67_z = as.numeric(extract_gene_z(expr_tcga_avg, "MKI67")),
  CellCycle_z = as.numeric(extract_gene_z(expr_tcga_avg, cell_cycle_genes))
)

tcga_adj <- tcga_risk %>%
  left_join(tcga_prolif, by = "id") %>%
  filter(!is.na(MKI67_z), !is.na(CellCycle_z), !is.na(time_months), !is.na(status), !is.na(risk_score))

cox_metric <- function(formula, data, label) {
  fit <- coxph(formula, data = data)
  tibble(
    Model = label,
    N = nrow(data),
    Events = sum(data$status == 1),
    C_index = as.numeric(summary(fit)$concordance[1]),
    C_index_SE = as.numeric(summary(fit)$concordance[2]),
    LogLik = as.numeric(logLik(fit)),
    AIC = AIC(fit)
  )
}

prolif_bench <- bind_rows(
  cox_metric(Surv(time_months, status) ~ risk_score, tcga_adj, "Nine-gene risk score"),
  cox_metric(Surv(time_months, status) ~ MKI67_z, tcga_adj, "MKI67 baseline"),
  cox_metric(Surv(time_months, status) ~ CellCycle_z, tcga_adj, "Cell-cycle metagene baseline"),
  cox_metric(Surv(time_months, status) ~ risk_score + MKI67_z, tcga_adj, "Risk score + MKI67"),
  cox_metric(Surv(time_months, status) ~ risk_score + CellCycle_z, tcga_adj, "Risk score + cell-cycle")
)

adj_fit <- coxph(Surv(time_months, status) ~ risk_score + MKI67_z + CellCycle_z, data = tcga_adj)
adj_sum <- summary(adj_fit)
prolif_effect <- tibble(
  Term = rownames(adj_sum$coefficients),
  HR = adj_sum$conf.int[, "exp(coef)"],
  CI95_L = adj_sum$conf.int[, "lower .95"],
  CI95_U = adj_sum$conf.int[, "upper .95"],
  P_value = adj_sum$coefficients[, "Pr(>|z|)"]
)

prolif_corr <- tibble(
  Metric = c("Risk score vs MKI67", "Risk score vs cell-cycle metagene", "MKI67 vs cell-cycle metagene"),
  Spearman_rho = c(
    suppressWarnings(cor(tcga_adj$risk_score, tcga_adj$MKI67_z, method = "spearman")),
    suppressWarnings(cor(tcga_adj$risk_score, tcga_adj$CellCycle_z, method = "spearman")),
    suppressWarnings(cor(tcga_adj$MKI67_z, tcga_adj$CellCycle_z, method = "spearman"))
  ),
  P_value = c(
    suppressWarnings(cor.test(tcga_adj$risk_score, tcga_adj$MKI67_z, method = "spearman")$p.value),
    suppressWarnings(cor.test(tcga_adj$risk_score, tcga_adj$CellCycle_z, method = "spearman")$p.value),
    suppressWarnings(cor.test(tcga_adj$MKI67_z, tcga_adj$CellCycle_z, method = "spearman")$p.value)
  )
)

write.csv(prolif_bench, file.path(res_dir, "tcga_proliferation_baseline_benchmark.csv"), row.names = FALSE)
write.csv(prolif_effect, file.path(res_dir, "tcga_risk_adjusted_for_proliferation.csv"), row.names = FALSE)
write.csv(prolif_corr, file.path(res_dir, "tcga_risk_proliferation_correlation.csv"), row.names = FALSE)

prediction_eval <- calibration_df %>%
  left_join(safe_read_csv(file.path(res_dir, "external_validation_stats.csv")) %>% dplyr::select(Dataset, C_index, AUC_3y), by = "Dataset") %>%
  transmute(
    dataset_id = Dataset,
    split_or_cohort = ifelse(Dataset == "GSE14520", "discovery", "external"),
    auc_3y = AUC_3y,
    c_index = C_index,
    calibration_intercept_3y = Calibration_intercept,
    calibration_slope_3y = Calibration_slope,
    brier_3y = Brier_3y_IPCW,
    n = N_total,
    n_landmark_3y = N_landmark
  )
write.csv(prediction_eval, file.path(bench_dir, "prediction_eval.csv"), row.names = FALSE)
write.table(
  prediction_eval,
  file.path(bench_dir, "prediction_eval.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

method_benchmark <- bind_rows(
  safe_read_csv(file.path(res_dir, "tcga_incremental_value.csv")) %>%
    transmute(Dataset, Model, N = Sample_N, Events, C_index, C_index_SE, LogLik, AIC),
  prolif_bench %>%
    mutate(Dataset = "TCGA-LIHC") %>%
    dplyr::select(Dataset, Model, N, Events, C_index, C_index_SE, LogLik, AIC)
) %>%
  distinct(Dataset, Model, .keep_all = TRUE)
write.csv(method_benchmark, file.path(bench_dir, "method_benchmark.csv"), row.names = FALSE)
write.table(
  method_benchmark,
  file.path(bench_dir, "method_benchmark.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

claim_effects <- bind_rows(
  safe_read_csv(file.path(res_dir, "risk_score_continuous_effect.csv")) %>%
    transmute(
      claim_id = "risk_score_survival",
      dataset_id = Dataset,
      outcome = "overall_survival",
      model = "cox_risk_score_per_1sd",
      effect_type = "hazard_ratio",
      effect = HR_per_1SD,
      ci_lower = CI95_L,
      ci_upper = CI95_U,
      pvalue = P_value,
      fdr = p.adjust(P_value, method = "BH"),
      n = Sample_N
    ),
  prolif_effect %>%
    filter(Term == "risk_score") %>%
    transmute(
      claim_id = "risk_score_adjusted_for_proliferation",
      dataset_id = "TCGA-LIHC",
      outcome = "overall_survival",
      model = "cox_risk_score_mki67_cellcycle",
      effect_type = "hazard_ratio",
      effect = HR,
      ci_lower = CI95_L,
      ci_upper = CI95_U,
      pvalue = P_value,
      fdr = P_value,
      n = nrow(tcga_adj)
    )
)
write.table(
  claim_effects,
  file.path(effect_dir, "claim_effects.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("[Model validation] Add BH-FDR columns to existing summary tables...")

add_fdr_if_possible <- function(infile, outfile, p_col_candidates) {
  df <- safe_read_csv(infile)
  p_col <- intersect(p_col_candidates, names(df))[1]
  if (!is.na(p_col)) {
    df$FDR_BH <- p.adjust(as.numeric(df[[p_col]]), method = "BH")
  }
  write.csv(df, outfile, row.names = FALSE)
}

add_fdr_if_possible(
  file.path(res_dir, "tcga_multiomics_keygene_mutation_assoc.csv"),
  file.path(res_dir, "tcga_multiomics_keygene_mutation_assoc_fdr.csv"),
  c("fisher_p", "P_value", "pvalue")
)
add_fdr_if_possible(
  file.path(res_dir, "immune_risk_group_diff.csv"),
  file.path(res_dir, "immune_risk_group_diff_fdr.csv"),
  c("P_value", "pvalue")
)
add_fdr_if_possible(
  file.path(res_dir, "checkpoint_risk_diff.csv"),
  file.path(res_dir, "checkpoint_risk_diff_fdr.csv"),
  c("P_value", "pvalue")
)

message("[Model validation] Additional validation and sensitivity analyses finished.")
