#!/usr/bin/env Rscript

# 02c_prognostic_model.R
# Ferroptosis-related prognostic model (DEG ∩ HCC-context ferroptosis genes)

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(glmnet)
  library(timeROC)
  library(cowplot)
  library(hgu133a2.db)
  library(AnnotationDbi)
})

proc_dir <- "data/processed"
ref_dir <- "data/references"
res_dir <- "results"
plot_dir <- "plots/supplementary"

dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

message("[预后模型] 开始构建铁死亡相关预后模型...")

expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
clinical_14 <- readRDS(file.path(proc_dir, "GSE14520_tumor_clinical.rds"))

deg_path <- file.path(res_dir, "deg_GSE14520_all.csv")
if (!file.exists(deg_path)) stop("缺少 deg_GSE14520_all.csv，请先运行 02b")

# Ferroptosis 基因集：主线使用 HCC-context（由 00c 生成），缺省回退到参考 FerrDb 列表
ferro_context_path <- file.path(res_dir, "ferroptosis_genes_hcc_context.csv")
ferro_ref_path <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")
ferro_path <- if (file.exists(ferro_context_path)) ferro_context_path else ferro_ref_path
if (!file.exists(ferro_path)) stop("缺少 ferroptosis 基因集：请先运行 00c 或准备 data/references/ferroptosis_genes_expanded.csv")

# DEG筛选（保持与主线一致）
deg_all <- read.csv(deg_path)
# DEG筛选（用于建模的候选池：严格）
deg_sig <- deg_all %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)

deg_genes <- toupper(deg_sig$Gene)
ferro_genes <- toupper(read.csv(ferro_path)$Gene)

# 识别景观基因 (DEG_strict ∩ Ferroptosis) 用于建模
landscape_genes <- intersect(unique(deg_genes), unique(ferro_genes))
message("[模型升级] 铁死亡景观基因数 (用于建模候选): ", length(landscape_genes))

if (length(landscape_genes) < 2) {
  stop("铁死亡差异景观基因为空，无法构建模型")
}

# 探针ID转基因符号
probe_ids <- rownames(expr_14)
gene_symbols <- mapIds(hgu133a2.db, keys = probe_ids,
                       column = "SYMBOL", keytype = "PROBEID",
                       multiVals = "first")

expr_gene <- expr_14
rownames(expr_gene) <- gene_symbols
expr_gene <- expr_gene[!is.na(rownames(expr_gene)), ]

# 合并重复基因
expr_df <- as.data.frame(expr_gene)
expr_df$gene <- rownames(expr_df)
expr_agg <- aggregate(. ~ gene, data = expr_df, FUN = mean)
rownames(expr_agg) <- expr_agg$gene
expr_agg$gene <- NULL
expr_gene <- as.matrix(expr_agg)

matched_gsm <- intersect(colnames(expr_gene), clinical_14$Affy_GSM)
expr_matched <- expr_gene[, matched_gsm]
clinical_matched <- clinical_14[match(matched_gsm, clinical_14$Affy_GSM), ]

zscore_by_gene <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

expr_scaled <- zscore_by_gene(expr_matched)

bootstrap_cindex <- function(df, n_boot = 500, seed = 123) {
  set.seed(seed)
  n <- nrow(df)
  vals <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    boot <- df[idx, , drop = FALSE]
    fit <- try(coxph(Surv(time, status) ~ risk_score, data = boot), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      vals[i] <- summary(fit)$concordance[1]
    }
  }
  vals <- vals[!is.na(vals)]
  ci <- quantile(vals, c(0.025, 0.975), na.rm = TRUE)
  list(
    mean = mean(vals),
    sd = sd(vals),
    ci_low = ci[[1]],
    ci_high = ci[[2]],
    n = length(vals)
  )
}

surv_data <- data.frame(
  sample = matched_gsm,
  time = as.numeric(clinical_matched$Survival.months),
  status = as.numeric(clinical_matched$Survival.status),
  stringsAsFactors = FALSE
)
surv_data <- surv_data[complete.cases(surv_data), ]

# 使用景观基因进行后续筛选
hub_in_expr <- intersect(landscape_genes, rownames(expr_scaled))
if (length(hub_in_expr) < 1) {
  stop("候选基因在表达矩阵中缺失")
}

# 单因素Cox筛选
cox_results <- data.frame(Gene = character(), HR = numeric(), HR_lower = numeric(),
                          HR_upper = numeric(), P_value = numeric(), stringsAsFactors = FALSE)

for (gene in hub_in_expr) {
  gene_expr <- as.numeric(expr_scaled[gene, surv_data$sample])
  if (sd(gene_expr, na.rm = TRUE) > 0.1) {
    tryCatch({
      cox_fit <- coxph(Surv(time, status) ~ gene_expr,
                       data = data.frame(time = surv_data$time,
                                        status = surv_data$status,
                                        gene_expr = gene_expr))
      summary_cox <- summary(cox_fit)
      cox_results <- rbind(cox_results, data.frame(
        Gene = gene,
        HR = summary_cox$conf.int[1, 1],
        HR_lower = summary_cox$conf.int[1, 3],
        HR_upper = summary_cox$conf.int[1, 4],
        P_value = summary_cox$coefficients[1, 5]
      ))
    }, error = function(e) {})
  }
}

cox_results <- cox_results %>% arrange(P_value)
write.csv(cox_results, file.path(res_dir, "cox_univariate.csv"), row.names = FALSE)

sig_genes <- cox_results %>% filter(P_value < 0.1) %>% pull(Gene)
if (length(sig_genes) < 2) {
  sig_genes <- head(cox_results$Gene, min(5, nrow(cox_results)))
}

if (length(sig_genes) < 1) {
  stop("Cox候选基因不足")
}

# LASSO-Cox
selected_genes <- sig_genes
if (length(sig_genes) >= 2) {
  x <- t(expr_scaled[sig_genes, surv_data$sample])
  y <- Surv(surv_data$time, surv_data$status)
  cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
  best_lambda <- cv_fit$lambda.min
  lasso_coef <- coef(cv_fit, s = best_lambda)
  lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)]
  if (length(lasso_genes) > 0) {
    selected_genes <- lasso_genes
  }

  # LASSO selection diagnostics (coefficient paths + CV deviance)
  # Saved for reviewer transparency (Supplementary Figure S8)
  pdf(file.path(plot_dir, "Supp_Figure_S8_LASSO_Selection.pdf"), width = 11, height = 5.5)
  par(mfrow = c(1, 2), mar = c(4.8, 4.4, 5.0, 1.2), mgp = c(2.6, 0.8, 0))
  plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "")
  title("LASSO-Cox coefficient paths", line = 3.4, cex.main = 0.95)
  abline(v = log(cv_fit$lambda.min), lty = 2, col = "#2F4F5F")
  plot(cv_fit, main = "")
  title("Cross-validated partial likelihood deviance", line = 3.4, cex.main = 0.95)
  dev.off()

  png(file.path(plot_dir, "Supp_Figure_S8_LASSO_Selection.png"), width = 3300, height = 1650, res = 300)
  par(mfrow = c(1, 2), mar = c(4.8, 4.4, 5.0, 1.2), mgp = c(2.6, 0.8, 0))
  plot(cv_fit$glmnet.fit, xvar = "lambda", label = FALSE, main = "")
  title("LASSO-Cox coefficient paths", line = 3.4, cex.main = 0.95)
  abline(v = log(cv_fit$lambda.min), lty = 2, col = "#2F4F5F")
  plot(cv_fit, main = "")
  title("Cross-validated partial likelihood deviance", line = 3.4, cex.main = 0.95)
  dev.off()
}

# 多因素 Cox 获取系数
if (length(selected_genes) > 1) {
  formula_str <- paste("Surv(time, status) ~", paste(paste0("`", selected_genes, "`"), collapse = " + "))
  model_data <- data.frame(surv_data, t(expr_scaled[selected_genes, surv_data$sample]))
  multi_cox <- coxph(as.formula(formula_str), data = model_data)
  coef_df <- data.frame(Gene = names(coef(multi_cox)), Coefficient = as.numeric(coef(multi_cox)))
} else {
  coef_df <- data.frame(Gene = selected_genes, Coefficient = 1)
}

write.csv(coef_df, file.path(res_dir, "prognostic_model_coef.csv"), row.names = FALSE)

# Risk score
risk_score <- rep(0, nrow(surv_data))
for (i in 1:nrow(coef_df)) {
  gene <- gsub("`", "", coef_df$Gene[i])
  coef <- coef_df$Coefficient[i]
  if (gene %in% rownames(expr_scaled)) {
    risk_score <- risk_score + coef * as.numeric(expr_scaled[gene, surv_data$sample])
  }
}

risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")

risk_df <- data.frame(
  sample = surv_data$sample,
  time = surv_data$time,
  status = surv_data$status,
  risk_score = risk_score,
  risk_group = factor(risk_group, levels = c("Low", "High"))
)

write.csv(risk_df, file.path(res_dir, "risk_score_data.csv"), row.names = FALSE)

# Model stats
surv_obj <- Surv(risk_df$time, risk_df$status)
cox_fit <- coxph(surv_obj ~ risk_score, data = risk_df)
cindex <- summary(cox_fit)$concordance[1]
logrank <- survdiff(surv_obj ~ risk_group, data = risk_df)
logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

boot_n <- suppressWarnings(as.integer(Sys.getenv("BOOTSTRAP_N", "500")))
if (is.na(boot_n) || boot_n < 50) boot_n <- 500
message("[预后模型] Bootstrap C-index (n=", boot_n, ")...")
boot_stats <- bootstrap_cindex(risk_df, n_boot = boot_n, seed = 123)
message(
  "[预后模型] Bootstrap C-index: ",
  sprintf("%.3f", boot_stats$mean),
  " (95% CI ",
  sprintf("%.3f", boot_stats$ci_low),
  "-",
  sprintf("%.3f", boot_stats$ci_high),
  ")"
)

roc_res <- timeROC(T = risk_df$time, delta = risk_df$status, marker = risk_df$risk_score,
                   cause = 1, weighting = "marginal", times = c(12, 36, 60), iid = FALSE)
auc_1 <- roc_res$AUC[roc_res$times == 12]
auc_3 <- roc_res$AUC[roc_res$times == 36]
auc_5 <- roc_res$AUC[roc_res$times == 60]

stats <- data.frame(
  Metric = c(
    "Log-rank p-value",
    "1-year AUC",
    "3-year AUC",
    "5-year AUC",
    "C-index",
    "Bootstrap C-index (mean)",
    "Bootstrap C-index (sd)",
    "Bootstrap C-index (95% CI lower)",
    "Bootstrap C-index (95% CI upper)",
    "Bootstrap N",
    "Selected Genes"
  ),
  Value = c(
    sprintf("%.3e", logrank_p),
    sprintf("%.3f", auc_1),
    sprintf("%.3f", auc_3),
    sprintf("%.3f", auc_5),
    sprintf("%.3f", cindex),
    sprintf("%.3f", boot_stats$mean),
    sprintf("%.3f", boot_stats$sd),
    sprintf("%.3f", boot_stats$ci_low),
    sprintf("%.3f", boot_stats$ci_high),
    sprintf("%d", boot_stats$n),
    paste(selected_genes, collapse = ", ")
  )
)
write.csv(stats, file.path(res_dir, "prognostic_model_stats.csv"), row.names = FALSE)

# Supplementary plot
km_fit <- survfit(Surv(time, status) ~ risk_group, data = risk_df)
p_km <- ggsurvplot(km_fit, data = risk_df, pval = TRUE, risk.table = FALSE,
                   palette = c("#2E9FDF", "#E7B800"),
                   xlab = "Time (months)", ylab = "Overall Survival",
                   legend.title = "", legend.labs = c("Low", "High"),
                   ggtheme = theme_bw(base_size = 10))$plot +
  labs(title = "Ferroptosis Signature: KM")

roc_df <- data.frame()
for (t in roc_res$times) {
  idx <- which(roc_res$times == t)
  tmp <- data.frame(
    FPR = roc_res$FP[, idx],
    TPR = roc_res$TP[, idx],
    Time = paste0(t / 12, "-year (AUC=", round(roc_res$AUC[idx], 3), ")")
  )
  roc_df <- rbind(roc_df, tmp)
}

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Time)) +
  geom_path(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Ferroptosis Signature: ROC", x = "1 - Specificity", y = "Sensitivity", color = "") +
  theme_bw(base_size = 10)

fig <- cowplot::plot_grid(p_km, p_roc, ncol = 2, labels = c("A", "B"))

ggsave(file.path(plot_dir, "FigureS_prognostic_model.pdf"), fig, width = 9, height = 4.5)

message("[预后模型] 完成")
message("  results/cox_univariate.csv")
message("  results/prognostic_model_coef.csv")
message("  results/risk_score_data.csv")
message("  results/prognostic_model_stats.csv")
message("  plots/supplementary/FigureS_prognostic_model.pdf")
