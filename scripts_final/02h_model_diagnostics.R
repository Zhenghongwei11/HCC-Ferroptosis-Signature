#!/usr/bin/env Rscript

# 02h_model_diagnostics.R
# Model diagnostics and sanity checks (no new data):
# - Cox PH assumption check for risk score (and optionally genes)
# - Random-gene assignment sanity check (random signatures with same coefficients)
# Outputs:
# - results/model_diagnostics_ph.csv
# - results/model_diagnostics_random_signature.csv
# - results/model_diagnostics_rmst.csv
# - plots/supplementary/FigureS_model_diagnostics_random_signature.pdf

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(ggplot2)
  library(hgu133a2.db)
  library(AnnotationDbi)
})

if (!dir.exists("data/processed")) {
  stop("Please run from project root")
}

proc_dir <- "data/processed"
res_dir <- "results"
plot_dir <- "plots/supplementary"

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

coef_df <- read.csv(file.path(res_dir, "prognostic_model_coef.csv"))
coef_vec <- setNames(coef_df$Coefficient, coef_df$Gene)
model_genes <- names(coef_vec)

gene_synonyms <- list(
  PKM = c("PKM", "PKM2")
)

harmonize_model_gene_names <- function(expr_mat, model_genes, cohort_label = "cohort") {
  for (gene in model_genes) {
    if (gene %in% rownames(expr_mat)) next
    aliases <- gene_synonyms[[gene]]
    if (is.null(aliases)) next
    alias_hit <- aliases[aliases %in% rownames(expr_mat)]
    alias_hit <- setdiff(alias_hit, gene)
    if (length(alias_hit) > 0) {
      expr_mat <- rbind(expr_mat, expr_mat[alias_hit[1], , drop = FALSE])
      rownames(expr_mat)[nrow(expr_mat)] <- gene
      message(sprintf("[Diagnostics] %s: using %s as synonym for %s", cohort_label, alias_hit[1], gene))
    }
  }
  expr_mat
}

maybe_log2p1 <- function(expr_mat) {
  q99 <- suppressWarnings(as.numeric(quantile(as.numeric(expr_mat), 0.99, na.rm = TRUE)))
  if (!is.na(q99) && q99 > 100) {
    return(log2(pmax(expr_mat, 0) + 1))
  }
  expr_mat
}

build_gse14520_model_expr <- function() {
  expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
  probe_ids <- rownames(expr_14)
  gene_symbols <- mapIds(
    hgu133a2.db,
    keys = probe_ids,
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

zscore_by_gene <- function(expr_mat) {
  z <- t(scale(t(expr_mat)))
  z[is.na(z)] <- 0
  z
}

calc_risk_score <- function(expr_mat, coef_vec) {
  expr_mat <- harmonize_model_gene_names(expr_mat, names(coef_vec))
  missing <- setdiff(names(coef_vec), rownames(expr_mat))
  if (length(missing) > 0) stop("Missing model genes: ", paste(missing, collapse = ", "))
  common <- names(coef_vec)
  colSums(expr_mat[common, , drop = FALSE] * coef_vec[common])
}

# Load discovery cohort (GSE14520) tumor samples using the same probe-to-gene
# mapping as model development.
train_expr <- build_gse14520_model_expr()
train_expr <- maybe_log2p1(train_expr)

# Use the already-curated survival table used by the main pipeline.
train_surv <- read.csv(file.path(res_dir, "risk_score_data.csv")) %>%
  transmute(id = sample, time_months = as.numeric(time), status = as.numeric(status)) %>%
  filter(!is.na(time_months), !is.na(status))

train_ids <- intersect(colnames(train_expr), train_surv$id)
train_expr <- train_expr[, train_ids, drop = FALSE]
train_surv <- train_surv %>% filter(id %in% train_ids)
train_surv <- train_surv[match(train_ids, train_surv$id), ]

# Load TCGA-LIHC
expr_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_expr.rds"))
clin_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_clinical.rds"))

# Average TCGA aliquots per patient
patient_id <- substr(colnames(expr_tcga), 1, 12)
expr_tcga_t <- t(expr_tcga)
expr_tcga_sum <- rowsum(expr_tcga_t, patient_id)
expr_tcga_avg <- expr_tcga_sum / as.vector(table(patient_id))
expr_tcga_avg <- t(expr_tcga_avg)

clin_tcga$patient_id <- clin_tcga$bcr_patient_barcode
clin_tcga$time_days <- ifelse(!is.na(clin_tcga$days_to_death), clin_tcga$days_to_death, clin_tcga$days_to_last_follow_up)
clin_tcga$status <- ifelse(tolower(clin_tcga$vital_status) == "dead", 1, 0)
clin_tcga$time_months <- as.numeric(clin_tcga$time_days) / 30.4

# Helper to build df for a cohort
build_df <- function(expr_mat, clin_df, id_col, time_col, status_col) {
  # z-score across genes within cohort on the intersected genes
  z <- zscore_by_gene(expr_mat)
  risk <- calc_risk_score(z, coef_vec)
  df <- data.frame(id = names(risk), risk_score = as.numeric(risk))
  clin_df2 <- clin_df[, c(id_col, time_col, status_col)]
  colnames(clin_df2) <- c("id", "time_months", "status")
  out <- merge(df, clin_df2, by = "id")
  out <- out %>% filter(!is.na(time_months), !is.na(status))
  out
}

train_df <- build_df(train_expr, train_surv, "id", "time_months", "status")
tcga_df <- build_df(expr_tcga_avg, clin_tcga, "patient_id", "time_months", "status")

# PH assumption check for risk score in both cohorts
ph_check <- function(df, cohort_name) {
  fit <- coxph(Surv(time_months, status) ~ risk_score, data = df)
  zph <- cox.zph(fit)
  # zph$table has row for risk_score and GLOBAL
  tab <- as.data.frame(zph$table)
  tab$term <- rownames(tab)
  tab$cohort <- cohort_name
  rownames(tab) <- NULL
  tab
}

ph_tab <- bind_rows(
  ph_check(train_df, "GSE14520"),
  ph_check(tcga_df, "TCGA-LIHC")
) %>%
  dplyr::select(cohort, term, chisq, df, p)

write.csv(ph_tab, file.path(res_dir, "model_diagnostics_ph.csv"), row.names = FALSE)

# RMST sensitivity analysis (median-split groups; robust to PH violations)
rmst_from_survfit <- function(sf, tau) {
  tt <- sort(unique(c(0, sf$time[sf$time <= tau], tau)))
  s <- summary(sf, times = tt)
  t <- s$time
  surv <- s$surv
  if (length(t) < 2) return(NA_real_)
  dt <- diff(t)
  # survival is piecewise-constant over [t_i, t_{i+1}) with value surv_i
  sum(dt * surv[-length(surv)])
}

rmst_compare <- function(df, cohort_name, tau_months = 60, n_boot = 500, seed = 20260129) {
  set.seed(seed)
  df <- df %>% filter(!is.na(time_months), !is.na(status), !is.na(risk_score))
  df$risk_group <- ifelse(df$risk_score >= median(df$risk_score), "High", "Low")
  tau <- min(tau_months, max(df$time_months, na.rm = TRUE))
  if (!is.finite(tau) || tau <= 0) {
    return(data.frame(
      cohort = cohort_name, tau = NA_real_, rmst_low = NA_real_, rmst_high = NA_real_,
      diff_high_minus_low = NA_real_, ci_low = NA_real_, ci_high = NA_real_
    ))
  }

  calc_once <- function(d) {
    sf_low <- survfit(Surv(time_months, status) ~ 1, data = d %>% filter(risk_group == "Low"))
    sf_high <- survfit(Surv(time_months, status) ~ 1, data = d %>% filter(risk_group == "High"))
    rl <- rmst_from_survfit(sf_low, tau)
    rh <- rmst_from_survfit(sf_high, tau)
    c(rmst_low = rl, rmst_high = rh, diff = rh - rl)
  }

  obs <- calc_once(df)
  diffs <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    boot <- df[idx, , drop = FALSE]
    if (length(unique(boot$risk_group)) < 2) next
    diffs[i] <- calc_once(boot)[["diff"]]
  }
  diffs <- diffs[is.finite(diffs)]
  ci <- if (length(diffs) >= 50) quantile(diffs, c(0.025, 0.975), na.rm = TRUE) else c(NA_real_, NA_real_)

  data.frame(
    cohort = cohort_name,
    tau = tau,
    rmst_low = unname(obs[["rmst_low"]]),
    rmst_high = unname(obs[["rmst_high"]]),
    diff_high_minus_low = unname(obs[["diff"]]),
    ci_low = as.numeric(ci[[1]]),
    ci_high = as.numeric(ci[[2]])
  )
}

rmst_tab <- bind_rows(
  rmst_compare(train_df, "GSE14520"),
  rmst_compare(tcga_df, "TCGA-LIHC")
)
write.csv(rmst_tab, file.path(res_dir, "model_diagnostics_rmst.csv"), row.names = FALSE)

# Random signature sanity check:
# Keep the same coefficient vector values, but assign them to random genes.
# This tests whether comparable performance could be achieved by arbitrary genes.
random_signature_test <- function(df_expr, df_surv, cohort_name, n_iter = 1000, seed = 20260129) {
  set.seed(seed)

  # Gene universe: genes with variance > 0
  v <- apply(df_expr, 1, var, na.rm = TRUE)
  universe <- names(v)[is.finite(v) & v > 0]

  k <- length(model_genes)
  coef_vals <- unname(coef_vec)
  names(coef_vals) <- model_genes

  # Observed
  z <- zscore_by_gene(df_expr)
  risk_obs <- calc_risk_score(z, coef_vec)
  df0 <- data.frame(risk_score = as.numeric(risk_obs), time_months = df_surv$time_months, status = df_surv$status)
  fit0 <- coxph(Surv(time_months, status) ~ risk_score, data = df0)
  c_obs <- summary(fit0)$concordance[1]

  c_rand <- rep(NA_real_, n_iter)
  for (i in seq_len(n_iter)) {
    genes_i <- sample(universe, k, replace = FALSE)
    coef_i <- coef_vals
    names(coef_i) <- genes_i
    common <- intersect(rownames(df_expr), names(coef_i))
    z_i <- zscore_by_gene(df_expr[common, , drop = FALSE])
    risk_i <- calc_risk_score(z_i, coef_i)
    df_i <- data.frame(risk_score = as.numeric(risk_i), time_months = df_surv$time_months, status = df_surv$status)
    fit_i <- coxph(Surv(time_months, status) ~ risk_score, data = df_i)
    c_rand[i] <- summary(fit_i)$concordance[1]
  }

  p_emp <- (sum(c_rand >= c_obs, na.rm = TRUE) + 1) / (sum(!is.na(c_rand)) + 1)

  data.frame(
    cohort = cohort_name,
    cindex_observed = as.numeric(c_obs),
    cindex_random_mean = mean(c_rand, na.rm = TRUE),
    cindex_random_median = median(c_rand, na.rm = TRUE),
    cindex_random_q05 = as.numeric(quantile(c_rand, 0.05, na.rm = TRUE)),
    cindex_random_q95 = as.numeric(quantile(c_rand, 0.95, na.rm = TRUE)),
    p_empirical_ge = p_emp
  )
}

# Align survival vectors to expression columns for both cohorts
train_surv <- train_df %>% arrange(match(id, colnames(train_expr)))
train_expr2 <- train_expr[, train_surv$id, drop = FALSE]
train_surv2 <- data.frame(time_months = train_surv$time_months, status = train_surv$status)

# tcga
# tcga_df$id are patient IDs; expr_tcga_avg cols are patient IDs
idx_tcga <- match(tcga_df$id, colnames(expr_tcga_avg))
expr_tcga2 <- expr_tcga_avg[, idx_tcga, drop = FALSE]
tcga_surv2 <- data.frame(time_months = tcga_df$time_months, status = tcga_df$status)

message("[Diagnostics] Random signature sanity check (this may take a few minutes)...")
rand_tab <- bind_rows(
  random_signature_test(train_expr2, train_surv2, "GSE14520", n_iter = 1000),
  random_signature_test(expr_tcga2, tcga_surv2, "TCGA-LIHC", n_iter = 1000)
)

write.csv(rand_tab, file.path(res_dir, "model_diagnostics_random_signature.csv"), row.names = FALSE)

# Plot distribution for each cohort
plot_rand_hist <- function(df_expr, df_surv, cohort_name, n_iter = 1000, seed = 20260129) {
  set.seed(seed)
  v <- apply(df_expr, 1, var, na.rm = TRUE)
  universe <- names(v)[is.finite(v) & v > 0]
  k <- length(model_genes)
  coef_vals <- unname(coef_vec)

  z <- zscore_by_gene(df_expr)
  risk_obs <- calc_risk_score(z, coef_vec)
  df0 <- data.frame(risk_score = as.numeric(risk_obs), time_months = df_surv$time_months, status = df_surv$status)
  fit0 <- coxph(Surv(time_months, status) ~ risk_score, data = df0)
  c_obs <- summary(fit0)$concordance[1]

  c_rand <- rep(NA_real_, n_iter)
  for (i in seq_len(n_iter)) {
    genes_i <- sample(universe, k, replace = FALSE)
    coef_i <- coef_vals
    names(coef_i) <- genes_i
    common <- intersect(rownames(df_expr), names(coef_i))
    z_i <- zscore_by_gene(df_expr[common, , drop = FALSE])
    risk_i <- calc_risk_score(z_i, coef_i)
    df_i <- data.frame(risk_score = as.numeric(risk_i), time_months = df_surv$time_months, status = df_surv$status)
    fit_i <- coxph(Surv(time_months, status) ~ risk_score, data = df_i)
    c_rand[i] <- summary(fit_i)$concordance[1]
  }

  tibble(cindex = c_rand, cohort = cohort_name, observed = as.numeric(c_obs))
}

h1 <- plot_rand_hist(train_expr2, train_surv2, "GSE14520", n_iter = 1000)
h2 <- plot_rand_hist(expr_tcga2, tcga_surv2, "TCGA-LIHC", n_iter = 1000)
h <- bind_rows(h1, h2)

p <- ggplot(h, aes(x = cindex)) +
  geom_histogram(bins = 40, fill = "#D9EAF7", color = "#1F4E79") +
  geom_vline(aes(xintercept = observed), color = "#B22222", linewidth = 1) +
  facet_wrap(~ cohort, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Random-gene assignment sanity check (C-index)",
    x = "C-index (random signatures with same coefficient magnitudes)",
    y = "Count"
  )

ggsave(file.path(plot_dir, "FigureS_model_diagnostics_random_signature.pdf"), p, width = 10, height = 4.5)
ggsave(file.path(plot_dir, "FigureS_model_diagnostics_random_signature.png"), p, width = 10, height = 4.5, dpi = 300, bg = "white")

message("[Diagnostics] Done")
message("  ✅ results/model_diagnostics_ph.csv")
message("  ✅ results/model_diagnostics_random_signature.csv")
message("  ✅ plots/supplementary/FigureS_model_diagnostics_random_signature.pdf")
message("  ✅ plots/supplementary/FigureS_model_diagnostics_random_signature.png")
