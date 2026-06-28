#!/usr/bin/env Rscript
# 07_final_publication_figures.R - Publication-quality figure panels (Figure 2–5)

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(rms)
  library(timeROC)
  library(gridExtra)
  library(grid)
  library(cowplot)
  library(ggrepel)
  library(ggsci)
})

# 目录设置
detect_base_dir <- function() {
  candidates <- c(".", "HCC_Ferroptosis_Project")
  for (d in candidates) {
    if (dir.exists(file.path(d, "results")) && dir.exists(file.path(d, "plots")) && dir.exists(file.path(d, "data/processed"))) return(d)
  }
  stop("无法定位项目根目录。")
}

base_dir <- detect_base_dir()
res_dir <- file.path(base_dir, "results")
proc_dir <- file.path(base_dir, "data/processed")
ref_dir <- file.path(base_dir, "data/references") 
pub_dir <- file.path(base_dir, "plots/publication")
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)
supp_dir <- file.path(base_dir, "plots/supplementary")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

# Publication theme
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "#2B3A3F"),
      axis.line = element_line(color = "#9FB1B6", linewidth = 0.4),
      panel.grid.major = element_line(color = "#E6ECEE"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

cool_palette <- c(
  "#2F4F5F", # slate teal
  "#6FA8B5", # muted teal
  "#8FB0B9", # blue gray
  "#AFC8D0", # pale teal
  "#C8D7DB"  # light gray teal
)

scale_color_cool <- function(...) {
  scale_color_manual(values = cool_palette, ...)
}

scale_fill_cool <- function(...) {
  scale_fill_manual(values = cool_palette, ...)
}

scale_fill_diverging <- function(...) {
  scale_fill_gradient2(low = "#2F4F5F", mid = "#F2F5F6", high = "#5F9EA6", ...)
}

format_immune_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\bB cell\\b", "B-cell") %>%
    str_replace_all("\\bCD4 T cell\\b", "CD4 T-cell") %>%
    str_replace_all("\\bCD8 T cell\\b", "CD8 T-cell") %>%
    str_replace_all("\\bMacrophage\\b", "macrophage") %>%
    str_replace_all("\\bMonocyte\\b", "monocyte") %>%
    str_replace_all("\\bType 17 T helper\\b", "Type 17 T-helper") %>%
    str_replace_all("\\bCentral memory\\b", "Central-memory")
}

lollipop_h <- function(df, x_col, y_col, color_col = NULL, title = "") {
  p <- ggplot(df, aes(x = !!sym(x_col), y = reorder(!!sym(y_col), !!sym(x_col)))) +
    geom_segment(aes(x = 0, xend = !!sym(x_col), y = !!sym(y_col), yend = !!sym(y_col)),
                 color = "#B0BCC0", linewidth = 0.5)
  if (is.null(color_col)) {
    p <- p + geom_point(color = "#2F4F5F", size = 2.6, alpha = 0.9, na.rm = TRUE)
  } else {
    p <- p + geom_point(aes(color = !!sym(color_col)), size = 2.6, alpha = 0.9, na.rm = TRUE)
  }
  p + theme_pub() + labs(title = title, x = NULL, y = NULL)
}

dotplot_h <- function(df, x_col, y_col, color_col = NULL, title = "") {
  p <- ggplot(df, aes(x = !!sym(x_col), y = reorder(!!sym(y_col), !!sym(x_col)))) +
    geom_vline(xintercept = 0, color = "#D5DDE0", linewidth = 0.4)
  if (is.null(color_col)) {
    p <- p + geom_point(color = "#2F4F5F", size = 3.0)
  } else {
    p <- p + geom_point(aes(color = !!sym(color_col)), size = 3.0)
  }
  p + theme_pub() + labs(title = title, x = NULL, y = NULL)
}

save_fig <- function(name, plot, w=12, h=10) {
  ggsave(file.path(pub_dir, paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(pub_dir, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  message(paste0("  ", name, " (PDF+PNG) saved."))
}

# ============================================================ 
# 加载数据
# ============================================================ 
risk_data <- read.csv(file.path(res_dir, "risk_score_data.csv"))
risk_data$risk_group <- factor(risk_data$risk_group, levels = c("Low", "High"))
logrank_fit <- survdiff(Surv(time, status) ~ risk_group, data = risk_data)
logrank_p <- formatC(1 - pchisq(logrank_fit$chisq, df = length(logrank_fit$n) - 1), format = "e", digits = 2)
model_coef <- read.csv(file.path(res_dir, "prognostic_model_coef.csv"))
model_stats <- read.csv(file.path(res_dir, "prognostic_model_stats.csv"))
uni_cox <- read.csv(file.path(res_dir, "univariate_cox_clinical.csv"))
multi_cox <- read.csv(file.path(res_dir, "multivariate_cox_clinical.csv"))
immune_corr <- read.csv(file.path(res_dir, "immune_risk_correlation.csv"))
immune_score_file <- file.path(res_dir, "marker_immune_signature_scores.csv")
if (!file.exists(immune_score_file)) {
  stop("缺少 marker_immune_signature_scores.csv，请先运行 03a_immune_infiltration.R")
}
immune_scores <- read.csv(immune_score_file, row.names = 1)
checkpoint_diff <- read.csv(file.path(res_dir, "checkpoint_risk_diff.csv"))
clinical_14 <- readRDS(file.path(proc_dir, "GSE14520_tumor_clinical.rds"))

# ----------------------------------------------------------------------------- 
# Figure 2: DEG & Intersection
# ----------------------------------------------------------------------------- 
message("[2/5] Figure 2...")
deg_all <- read.csv(file.path(res_dir, "deg_GSE14520_all.csv"))
deg_all$Sig <- "NS"; deg_all$Sig[deg_all$adj.P.Val < 0.05 & deg_all$logFC > 1] <- "Up"; deg_all$Sig[deg_all$adj.P.Val < 0.05 & deg_all$logFC < -1] <- "Down"

p2a <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val+1e-300), color = Sig)) +
  geom_point(alpha = 0.3, size = 0.7) +
  scale_color_manual(values = c(Up = "#6FA8B5", Down = "#2F4F5F", NS = "#C8D7DB")) +
  geom_text_repel(data = head(deg_all %>% arrange(adj.P.Val), 10), aes(label = Gene), size = 3, fontface="bold", color="black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  theme_pub() + labs(title = "Volcano Plot (Tumor vs Non-tumor)", x = "log2FC", y = "-log10(FDR)")

expr_mat <- readRDS(file.path(proc_dir, "GSE14520_expr_symbol.rds"))
clinical_all <- readRDS(file.path(proc_dir, "GSE14520_clinical.rds"))
top_genes <- deg_all %>% arrange(adj.P.Val) %>% head(60) %>% pull(Gene)
top_genes <- intersect(top_genes, rownames(expr_mat)) %>% head(20)
expr_sub <- expr_mat[top_genes, , drop = FALSE]

group_info <- clinical_all %>%
  dplyr::select(geo_accession, `Tissue:ch1`) %>%
  mutate(Group = ifelse(grepl("Tumor", `Tissue:ch1`, ignore.case = TRUE), "Tumor", "Non-Tumor")) %>%
  dplyr::select(geo_accession, Group)

sample_order <- group_info$geo_accession[group_info$geo_accession %in% colnames(expr_sub)]
expr_long <- as.data.frame(expr_sub[, sample_order, drop = FALSE]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  group_by(Gene) %>%
  mutate(Z = as.numeric(scale(Expression))) %>%
  ungroup() %>%
  left_join(group_info, by = c("Sample" = "geo_accession"))

sample_order <- group_info %>%
  filter(geo_accession %in% sample_order) %>%
  arrange(Group) %>%
  pull(geo_accession)

expr_long$Sample <- factor(expr_long$Sample, levels = sample_order)
expr_long$Gene <- factor(expr_long$Gene, levels = rev(top_genes))

p2b <- ggplot(expr_long, aes(x = Sample, y = Gene, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Top DEGs Heatmap", x = NULL, y = NULL, fill = "Z-score")

ferro_context_path <- file.path(res_dir, "ferroptosis_genes_hcc_context.csv")
ferro_ref_path <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")
ferro_path <- if (file.exists(ferro_context_path)) ferro_context_path else ferro_ref_path
if (!file.exists(ferro_path)) stop("缺少 ferroptosis 基因集：请先运行 00c 或准备 data/references/ferroptosis_genes_expanded.csv")
ferro_genes <- read.csv(ferro_path)$Gene
deg_sig <- deg_all %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% pull(Gene)

deg_set <- unique(deg_sig)
ferro_set <- unique(ferro_genes)
deg_ferro <- intersect(deg_set, ferro_set)
only_deg <- setdiff(deg_set, ferro_set)

circle_points <- function(x0, y0, r, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(t), y = y0 + r * sin(t))
}

deg_circle <- circle_points(-1, 0, 1.6)
ferro_circle <- circle_points(1, 0, 1.6)

p2c <- ggplot() +
  geom_polygon(data = deg_circle, aes(x = x, y = y), fill = "#3C5488FF", alpha = 0.18, color = "#3C5488FF") +
  geom_polygon(data = ferro_circle, aes(x = x, y = y), fill = "#E64B35FF", alpha = 0.18, color = "#E64B35FF") +
  annotate("text", x = -2.2, y = -1.1, label = paste0("DEG\n(n=", length(deg_set), ")"), size = 3.6, color = "#3C5488FF") +
  annotate("text", x = 2.2, y = -1.1, label = paste0("Ferroptosis\n(n=", length(ferro_set), ")"), size = 3.6, color = "#E64B35FF") +
  annotate("text", x = -1.7, y = 0.2, label = length(only_deg), size = 4) +
  annotate("text", x = 1.7, y = 0.2, label = length(setdiff(ferro_set, deg_set)), size = 4) +
  annotate("text", x = 0, y = 0.2, label = length(deg_ferro), size = 5, fontface = "bold") +
  coord_equal() +
  theme_void() +
  labs(title = "HCC-specific ferroptosis candidates")

candidate_dat <- deg_all %>%
  filter(Gene %in% deg_ferro) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(logFC = mean(logFC, na.rm = TRUE), adj.P.Val = min(adj.P.Val, na.rm = TRUE), .groups = "drop") %>%
  arrange(adj.P.Val) %>%
  head(12) %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down"))
p2d <- lollipop_h(candidate_dat, "logFC", "Gene", "Direction", title = "Top HCC-context ferroptosis candidates") +
  scale_color_manual(values = c(Up = "#6FA8B5", Down = "#2F4F5F")) +
  labs(x = "log2FC", y = "") +
  theme(legend.position = "none")

fig2 <- plot_grid(p2a, p2b, p2c, p2d, labels = "AUTO", ncol = 2)
save_fig("Figure2_deg_ferroptosis", fig2)

# ----------------------------------------------------------------------------- 
# Figure 3: Prognostic Model (还原 A-D)
# ----------------------------------------------------------------------------- 
message("[3/5] Figure 3...")
km_fit <- survfit(Surv(time, status) ~ risk_group, data = risk_data)
km_obj <- ggsurvplot(
  km_fit,
  data = risk_data,
  pval = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  break.time.by = 12,
  palette = "npg",
  ggtheme = theme_pub(),
  legend.labs = c("Low Risk", "High Risk"),
  xlab = "Time (months)"
)
p3a_plot <- km_obj$plot +
  annotate("text", x = max(risk_data$time, na.rm = TRUE) * 0.6, y = 0.1, label = paste0("p = ", logrank_p), size = 3.2)
p3a_table <- km_obj$table + theme_pub()
p3a <- plot_grid(p3a_plot, p3a_table, ncol = 1, rel_heights = c(3, 1))

roc_times <- c(12, 36, 60)
roc_labels <- c("1-year", "3-year", "5-year")
roc_res <- timeROC(T = risk_data$time, delta = risk_data$status, marker = risk_data$risk_score, cause = 1, weighting = "marginal", times = roc_times, iid = FALSE)
roc_df <- data.frame()
for (i in seq_along(roc_times)) roc_df <- rbind(roc_df, data.frame(FPR=roc_res$FP[,i], TPR=roc_res$TP[,i], Time=roc_labels[i]))
p3b <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Time)) + geom_path(linewidth = 1.2) + scale_color_cool() + geom_abline(linetype="dashed") + theme_pub() + labs(title = "ROC Analysis")

risk_sorted <- risk_data %>% arrange(risk_score) %>% mutate(rank = row_number())
p3c <- ggplot(risk_sorted, aes(x = rank, y = risk_score, color = risk_group)) +
  geom_point(size = 1) +
  scale_color_cool(labels = c("Low risk", "High risk")) +
  theme_pub() +
  labs(title = "Risk Score Distribution", x = "Ranked patients", y = "Risk score", color = "Risk group")

model_coef$CoefficientSign <- ifelse(model_coef$Coefficient > 0, "Positive", "Negative")
p3d <- ggplot(model_coef, aes(x = Coefficient, y = reorder(Gene, Coefficient), fill = CoefficientSign)) + geom_col() + scale_fill_cool() + theme_pub() + labs(title = "Model Coefficients", y="", fill = "Coefficient sign")

fig3 <- plot_grid(p3a, p3b, p3c, p3d, labels = "AUTO", ncol = 2)
save_fig("Figure3_prognostic_model", fig3)

# ----------------------------------------------------------------------------- 
# Figure 4: Clinical Utility (还原 A-D)
# ----------------------------------------------------------------------------- 
message("[4/5] Figure 4...")
uni_cox <- uni_cox %>% mutate(Display = Variable)
multi_cox <- multi_cox %>%
  mutate(Display = case_when(
    Variable == "risk_score" ~ "Risk score",
    Variable == "tumor_sizeLarge" ~ "Tumor size (Large vs Small)",
    Variable == "afpHigh" ~ "AFP (High vs Low)",
    Variable == "tnm_groupAdvanced" ~ "TNM stage (Advanced vs Early)",
    TRUE ~ Variable
  ))

p4a <- ggplot(uni_cox, aes(x = HR, y = reorder(Display, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
  geom_point(aes(color = P_value < 0.05), size = 3) +
  scale_color_cool() + scale_x_log10() + theme_pub() +
  labs(title = "Univariate Cox", y = "Variables", x = "Hazard Ratio (log scale)")

p4b <- ggplot(multi_cox, aes(x = HR, y = reorder(Display, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#9FB1B6", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
  geom_point(aes(color = P_value < 0.05), size = 3) +
  scale_color_cool() + scale_x_log10() + theme_pub() +
  labs(title = "Multivariate Cox", y = "Variables", x = "Hazard Ratio (log scale)")

# DCA (Net benefit vs threshold probability)
# Landmark-style 3-year DCA:
# - Predicted risk: 3-year death probability from a Cox model using baseline cumulative hazard
# - Outcome: death within 36 months; samples censored before 36 months are excluded
full_data <- data.frame(time = risk_data$time, status = risk_data$status, risk = risk_data$risk_score)
cox_model <- coxph(Surv(time, status) ~ risk, data = full_data)
lp <- predict(cox_model, type = "lp")

t0 <- 36 # months
bh <- basehaz(cox_model, centered = FALSE)
H0_t0 <- approx(bh$time, bh$hazard, xout = t0, method = "linear", rule = 2)$y
pred_prob <- 1 - exp(-H0_t0 * exp(lp)) # P(death by 36 months)

eligible <- (full_data$time >= t0) | (full_data$status == 1 & full_data$time < t0)
dca_outcome <- ifelse(full_data$status == 1 & full_data$time <= t0, 1, 0)

pred_prob <- pred_prob[eligible]
dca_outcome <- dca_outcome[eligible]

thresholds <- seq(0.01, 0.99, by = 0.01)
calc_nb <- function(th) {
  pred_pos <- pred_prob >= th
  tp <- sum(pred_pos & dca_outcome == 1, na.rm = TRUE)
  fp <- sum(pred_pos & dca_outcome == 0, na.rm = TRUE)
  n <- sum(!is.na(pred_pos))
  (tp / n) - (fp / n) * (th / (1 - th))
}
nb_model <- sapply(thresholds, calc_nb)
event_rate <- mean(dca_outcome, na.rm = TRUE)
nb_all <- event_rate - (1 - event_rate) * (thresholds / (1 - thresholds))
nb_none <- rep(0, length(thresholds))

dca_df <- data.frame(
  Threshold = thresholds,
  Model = nb_model,
  TreatAll = nb_all,
  TreatNone = nb_none
)

dca_long <- dca_df %>%
  pivot_longer(cols = c(Model, TreatAll, TreatNone), names_to = "Strategy", values_to = "NetBenefit")

p4c_y <- dca_long %>%
  filter(Threshold <= 0.6) %>%
  summarise(ymin = min(NetBenefit, na.rm = TRUE), ymax = max(NetBenefit, na.rm = TRUE))
p4c <- ggplot(dca_long, aes(x = Threshold, y = NetBenefit, color = Strategy)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c(Model = "#2F4F5F", TreatAll = "#9FB1B6", TreatNone = "#C8D7DB")) +
  theme_pub() +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(p4c_y$ymin, p4c_y$ymax)) +
  labs(title = "Decision Curve Analysis (3-year)", x = "Threshold probability", y = "Net benefit")

calc_dca_long <- function(risk_data, t0_months, time_label) {
  full_data <- data.frame(time = risk_data$time, status = risk_data$status, risk = risk_data$risk_score)
  cox_model <- coxph(Surv(time, status) ~ risk, data = full_data)
  lp <- predict(cox_model, type = "lp")

  bh <- basehaz(cox_model, centered = FALSE)
  H0_t0 <- approx(bh$time, bh$hazard, xout = t0_months, method = "linear", rule = 2)$y
  pred_prob <- 1 - exp(-H0_t0 * exp(lp)) # P(death by t0)

  eligible <- (full_data$time >= t0_months) | (full_data$status == 1 & full_data$time < t0_months)
  dca_outcome <- ifelse(full_data$status == 1 & full_data$time <= t0_months, 1, 0)

  pred_prob <- pred_prob[eligible]
  dca_outcome <- dca_outcome[eligible]

  thresholds <- seq(0.01, 0.99, by = 0.01)
  calc_nb <- function(th) {
    pred_pos <- pred_prob >= th
    tp <- sum(pred_pos & dca_outcome == 1, na.rm = TRUE)
    fp <- sum(pred_pos & dca_outcome == 0, na.rm = TRUE)
    n <- sum(!is.na(pred_pos))
    (tp / n) - (fp / n) * (th / (1 - th))
  }
  nb_model <- sapply(thresholds, calc_nb)
  event_rate <- mean(dca_outcome, na.rm = TRUE)
  nb_all <- event_rate - (1 - event_rate) * (thresholds / (1 - thresholds))
  nb_none <- rep(0, length(thresholds))

  dca_df <- data.frame(
    Time = time_label,
    Threshold = thresholds,
    Model = nb_model,
    TreatAll = nb_all,
    TreatNone = nb_none
  )
  dca_df %>%
    pivot_longer(cols = c(Model, TreatAll, TreatNone), names_to = "Strategy", values_to = "NetBenefit")
}

dca_1y <- calc_dca_long(risk_data, t0_months = 12, time_label = "1-year")
dca_3y <- calc_dca_long(risk_data, t0_months = 36, time_label = "3-year")
dca_5y <- calc_dca_long(risk_data, t0_months = 60, time_label = "5-year")
dca_all_long <- bind_rows(dca_1y, dca_3y, dca_5y)
dca_all_plot <- dca_all_long %>% filter(Threshold <= 0.6)
dca_y <- dca_all_plot %>%
  summarise(ymin = min(NetBenefit, na.rm = TRUE), ymax = max(NetBenefit, na.rm = TRUE))

p_dca_multi <- ggplot(dca_all_plot, aes(x = Threshold, y = NetBenefit, color = Strategy)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c(Model = "#2F4F5F", TreatAll = "#9FB1B6", TreatNone = "#C8D7DB")) +
  theme_pub() +
  facet_wrap(~Time, ncol = 3) +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(dca_y$ymin, dca_y$ymax)) +
  labs(title = "Decision Curve Analysis (1/3/5-year)", x = "Threshold probability", y = "Net benefit")

ggsave(file.path(supp_dir, "Supp_Figure_S10_DCA_1_3_5y.pdf"), p_dca_multi, width = 13.5, height = 4.5, device = cairo_pdf)
ggsave(file.path(supp_dir, "Supp_Figure_S10_DCA_1_3_5y.png"), p_dca_multi, width = 13.5, height = 4.5, dpi = 300, bg = "white")
message("  Supplementary DCA (1/3/5-year) saved.")

nomo_file <- file.path(res_dir, "nomogram_model_performance.csv")
risk_cindex <- as.numeric(model_stats$Value[model_stats$Metric == "C-index"])
risk_ci_lower <- as.numeric(model_stats$Value[model_stats$Metric == "Bootstrap C-index (95% CI lower)"])
risk_ci_upper <- as.numeric(model_stats$Value[model_stats$Metric == "Bootstrap C-index (95% CI upper)"])
nomogram_cindex <- NA
nomo_ci_lower <- NA
nomo_ci_upper <- NA
if (file.exists(nomo_file)) {
  nomo <- read.csv(nomo_file)
  nomogram_cindex <- as.numeric(nomo$Value[nomo$Metric == "C-index"])
  nomo_se <- as.numeric(nomo$Value[nomo$Metric == "C-index SE"])
  if (!is.na(nomo_se)) {
    nomo_ci_lower <- max(0, nomogram_cindex - 1.96 * nomo_se)
    nomo_ci_upper <- min(1, nomogram_cindex + 1.96 * nomo_se)
  }
}
cindex_df <- data.frame(
  Model = c("Risk Score", "Nomogram"),
  Cindex = c(risk_cindex, nomogram_cindex),
  Lower = c(risk_ci_lower, nomo_ci_lower),
  Upper = c(risk_ci_upper, nomo_ci_upper)
)
cindex_df <- cindex_df[!is.na(cindex_df$Cindex), ]

p4d <- ggplot(cindex_df, aes(x = Cindex, y = Model)) +
  annotate("rect", xmin = 0.5, xmax = 0.8, ymin = -Inf, ymax = Inf, fill = "#F2F5F6") +
  geom_vline(xintercept = 0.5, color = "#9FB1B6", linewidth = 0.5, linetype = "dashed") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.15, linewidth = 0.6, color = "#6FA8B5", na.rm = TRUE) +
  geom_point(shape = 21, fill = "white", color = "#2F4F5F", size = 3.6) +
  geom_text(aes(label = sprintf("%.3f", Cindex)), hjust = -0.3, size = 3.2) +
  theme_pub() +
  labs(title = "C-index Comparison", x = "C-index", y = NULL) +
  theme(legend.position = "none")

fig4 <- plot_grid(p4a, p4b, p4c, p4d, labels = "AUTO", ncol = 2)
save_fig("Figure4_decision_curve_nomogram", fig4)

# ----------------------------------------------------------------------------- 
# Figure 5: Immune (还原 A-D)
# ----------------------------------------------------------------------------- 
message("[5/5] Figure 5...")
immune_corr_plot <- immune_corr %>%
  arrange(desc(abs(Correlation))) %>%
  head(12) %>%
  mutate(
    Direction = ifelse(Correlation >= 0, "Positive", "Negative"),
    Cell_Label = format_immune_label(Cell_Type)
  )
p5a <- lollipop_h(immune_corr_plot, "Correlation", "Cell_Label", "Direction", title = "Immune Correlation") +
  scale_color_manual(values = c(Positive = "#6FA8B5", Negative = "#2F4F5F")) +
  labs(x = "Correlation", y = NULL)

imm_box <- immune_scores %>%
  rownames_to_column("sample") %>%
  inner_join(risk_data, by="sample") %>%
  dplyr::select(risk_group, 1:5) %>%
  dplyr::select(-sample) %>%
  pivot_longer(-risk_group) %>%
  mutate(name = format_immune_label(name))
p5b <- ggplot(imm_box, aes(x=name, y=value, fill=risk_group)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_cool(labels = c("Low risk", "High risk")) +
  theme_pub() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune Infiltration", x = "Cell type", y = "Signature score", fill = "Risk group")

check_plot <- checkpoint_diff %>%
  arrange(P_value) %>%
  head(10) %>%
  mutate(Direction = ifelse(Log2FC >= 0, "Up in High", "Down in High"))
p5c <- ggplot(check_plot, aes(x = Log2FC, y = reorder(Gene, Log2FC), color = Direction, size = -log10(P_value))) +
  geom_vline(xintercept = 0, color = "#9FB1B6", linewidth = 0.5) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("Up in High" = "#6FA8B5", "Down in High" = "#2F4F5F")) +
  scale_size(range = c(2.5, 6)) +
  theme_pub() +
  labs(title = "Checkpoints (High vs Low)", x = "Log2FC (High - Low)", y = "Checkpoint", size = "-log10(P)")

merged_imm <- inner_join(risk_data, immune_scores %>% rownames_to_column("sample"), by="sample")
p5d <- ggplot(merged_imm, aes(x = risk_score, y = M2_Macrophage)) +
  geom_point(aes(color = risk_group), alpha = 0.5) +
  geom_smooth(method = "lm", color="black") +
  stat_cor(method="spearman") +
  scale_color_cool(labels = c("Low risk", "High risk")) +
  theme_pub() +
  labs(
    title = "Risk vs M2 macrophages",
    x = "Risk score",
    y = "M2 macrophage signature score",
    color = "Risk group"
  )

fig5 <- plot_grid(p5a, p5b, p5c, p5d, labels = "AUTO", ncol = 2)
save_fig("Figure5_immune_analysis", fig5)

message("\n============================================================")
message("所有主图已生成。输出目录: plots/publication/")
message("Main publication figures generated.")
message("============================================================")
