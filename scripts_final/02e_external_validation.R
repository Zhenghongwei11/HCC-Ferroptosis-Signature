#!/usr/bin/env Rscript

# 02e_external_validation.R
# 外部验证：TCGA-LIHC + GSE76427 + GSE10143 + (可选) GSE27150 + (可选) ICGC-LIRI-JP(HCCDB18)
# 输出：risk score、外部验证统计表、可选图版

# 禁用交互式图形设备
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GEOquery)
  library(survival)
  library(survminer)
  library(timeROC)
  library(cowplot)
  library(hgu133a2.db)
  library(AnnotationDbi)
})

if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
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
      message(sprintf("[外部验证] %s: using %s as synonym for %s", cohort_label, alias_hit[1], gene))
    }
  }
  expr_mat
}

calc_risk_score <- function(expr_mat, coef_vec, cohort_label = "cohort", strict = FALSE, min_genes = 2) {
  expr_mat <- harmonize_model_gene_names(expr_mat, names(coef_vec), cohort_label = cohort_label)
  missing <- setdiff(names(coef_vec), rownames(expr_mat))
  if (strict && length(missing) > 0) {
    stop("外部队列缺少模型基因，无法计算完整九基因风险评分: ", paste(missing, collapse = ", "))
  }
  common <- intersect(names(coef_vec), rownames(expr_mat))
  if (length(common) < min_genes) {
    stop("外部队列可用模型基因过少，无法计算风险评分: ", cohort_label)
  }
  if (length(missing) > 0) {
    message(sprintf("[外部验证] %s: using %d/%d model genes; missing %s",
                    cohort_label, length(common), length(coef_vec), paste(missing, collapse = ", ")))
  }
  score <- colSums(expr_mat[common, , drop = FALSE] * coef_vec[common])
  attr(score, "genes_used_n") <- length(common)
  attr(score, "genes_total_n") <- length(coef_vec)
  attr(score, "genes_used") <- paste(common, collapse = ";")
  attr(score, "missing_genes") <- if (length(missing) == 0) "" else paste(missing, collapse = ";")
  score
}

add_gene_coverage <- function(df, risk_score) {
  df$Genes_Used_N <- attr(risk_score, "genes_used_n")
  df$Genes_Total_N <- attr(risk_score, "genes_total_n")
  df$Genes_Used <- attr(risk_score, "genes_used")
  df$Missing_Model_Genes <- attr(risk_score, "missing_genes")
  df
}

map_expression_to_symbols <- function(es) {
  expr_raw <- exprs(es)
  feat <- fData(es)
  gene_col <- grep("gene symbol|gene_symbol|symbol", colnames(feat), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gene_col)) {
    stop("找不到基因符号列")
  }
  symbols <- feat[[gene_col]]
  symbols <- trimws(symbols)
  symbols <- sub(" ///.*$", "", symbols)
  symbols <- sub(" //.*$", "", symbols)
  symbols <- sub(";.*$", "", symbols)

  keep <- !is.na(symbols) & symbols != ""
  expr_raw <- expr_raw[keep, , drop = FALSE]
  symbols <- symbols[keep]

  expr_sum <- rowsum(expr_raw, symbols)
  expr_cnt <- rowsum(matrix(1, nrow(expr_raw), ncol(expr_raw)), symbols)
  expr_sum / expr_cnt
}

map_expression_with_gpl <- function(expr_raw, gpl_id) {
  gpl <- getGEO(gpl_id)
  anno <- Table(gpl)
  id_col <- grep("^ID$|^ID_REF$", colnames(anno), ignore.case = TRUE, value = TRUE)[1]
  gene_col <- grep("gene symbol|gene_symbol|symbol", colnames(anno), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col) || is.na(gene_col)) {
    stop("GPL 注释表缺少探针或基因符号列")
  }

  anno <- anno[, c(id_col, gene_col)]
  colnames(anno) <- c("Probe", "Symbol")
  anno$Symbol <- trimws(anno$Symbol)
  anno$Symbol <- sub(" ///.*$", "", anno$Symbol)
  anno$Symbol <- sub(" //.*$", "", anno$Symbol)
  anno$Symbol <- sub(";.*$", "", anno$Symbol)

  probe <- rownames(expr_raw)
  matched <- anno[match(probe, anno$Probe), ]
  symbols <- matched$Symbol

  keep <- !is.na(symbols) & symbols != ""
  expr_raw <- expr_raw[keep, , drop = FALSE]
  symbols <- symbols[keep]

  expr_sum <- rowsum(expr_raw, symbols)
  expr_cnt <- rowsum(matrix(1, nrow(expr_raw), ncol(expr_raw)), symbols)
  expr_sum / expr_cnt
}

zscore_by_gene <- function(expr_mat) {
  z <- t(scale(t(expr_mat)))
  z[is.na(z)] <- 0
  z
}

maybe_log2p1 <- function(expr_mat) {
  q99 <- suppressWarnings(as.numeric(quantile(as.numeric(expr_mat), 0.99, na.rm = TRUE)))
  if (!is.na(q99) && q99 > 100) {
    # guard against negative values from normalized matrices
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

train_frozen_scaler <- function(model_genes) {
  train_ids_path <- file.path(res_dir, "risk_score_data.csv")
  if (!file.exists(file.path(proc_dir, "GSE14520_expr.rds")) || !file.exists(train_ids_path)) {
    stop("缺少训练集表达或风险评分文件，无法构建冻结尺度参数")
  }
  train_expr <- build_gse14520_model_expr()
  train_expr <- maybe_log2p1(train_expr)
  train_ids <- read.csv(train_ids_path)$sample
  train_ids <- intersect(colnames(train_expr), train_ids)
  if (length(train_ids) < 20) stop("训练集样本过少，无法估计冻结尺度参数")
  train_expr <- harmonize_model_gene_names(train_expr, model_genes, cohort_label = "GSE14520")
  missing <- setdiff(model_genes, rownames(train_expr))
  if (length(missing) > 0) stop("训练集缺少模型基因，无法构建冻结尺度参数: ", paste(missing, collapse = ", "))
  genes <- model_genes
  mu <- rowMeans(train_expr[genes, train_ids, drop = FALSE], na.rm = TRUE)
  sdv <- apply(train_expr[genes, train_ids, drop = FALSE], 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(mu = mu, sd = sdv, genes = genes)
}

zscore_by_train <- function(expr_mat, scaler) {
  expr_mat <- maybe_log2p1(expr_mat)
  expr_mat <- harmonize_model_gene_names(expr_mat, scaler$genes)
  missing <- setdiff(scaler$genes, rownames(expr_mat))
  if (length(missing) > 0) {
    message("[外部验证] frozen scaling: missing model genes after synonym mapping: ", paste(missing, collapse = ", "))
  }
  common <- intersect(scaler$genes, rownames(expr_mat))
  if (length(common) < 2) stop("外部队列可用模型基因过少，无法冻结尺度标准化")
  z <- expr_mat[common, , drop = FALSE]
  z <- (z - scaler$mu[common]) / scaler$sd[common]
  z[is.na(z)] <- 0
  z
}

extract_characteristic <- function(pdat, pattern) {
  char_cols <- grep("characteristics", colnames(pdat), ignore.case = TRUE, value = TRUE)
  if (length(char_cols) == 0) return(rep(NA_character_, nrow(pdat)))
  pdat[char_cols] <- lapply(pdat[char_cols], as.character)
  apply(pdat[, char_cols, drop = FALSE], 1, function(x) {
    hit <- x[grepl(pattern, x, ignore.case = TRUE)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  })
}

extract_last_number <- function(x) {
  if (is.na(x)) return(NA_real_)
  m <- sub(".*:\\s*([0-9]+)\\s*$", "\\1", x)
  if (identical(m, x)) return(NA_real_)
  suppressWarnings(as.numeric(m))
}

extract_yes_no <- function(x) {
  if (is.na(x)) return(NA_character_)
  m <- sub(".*:\\s*([Yy]es|[Nn]o)\\b.*", "\\1", x)
  if (identical(m, x)) return(NA_character_)
  ifelse(tolower(m) == "yes", "Yes", "No")
}

collapse_characteristics <- function(pdat) {
  char_cols <- grep("^characteristics", colnames(pdat), ignore.case = TRUE, value = TRUE)
  if (length(char_cols) == 0) {
    return(rep("", nrow(pdat)))
  }
  apply(pdat[, char_cols, drop = FALSE], 1, function(x) {
    paste(na.omit(as.character(x)), collapse = ";")
  })
}

prepare_hccdb18_lirijp <- function(raw_dir = file.path("data", "raw", "HCCDB18")) {
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

  expr_txt <- file.path(raw_dir, "HCCDB18_mRNA_level3.txt")
  samp_txt <- file.path(raw_dir, "HCCDB18.sample.txt")
  pat_txt <- file.path(raw_dir, "HCCDB18.patient.txt")

  if (!file.exists(expr_txt) || !file.exists(samp_txt) || !file.exists(pat_txt)) {
    message("[HCCDB18] 检测到原始文件缺失，开始从 HCCDB 下载并解压...")
    expr_zip <- file.path(raw_dir, "HCCDB18_mRNA_level3.zip")
    samp_zip <- file.path(raw_dir, "HCCDB18.sample.zip")
    pat_zip <- file.path(raw_dir, "HCCDB18.patient.zip")

    download.file("http://lifeome.net/database/hccdb/download/HCCDB18_mRNA_level3.zip", expr_zip, mode = "wb", quiet = TRUE)
    download.file("http://lifeome.net/database/hccdb/download/sample/HCCDB18.sample.zip", samp_zip, mode = "wb", quiet = TRUE)
    download.file("http://lifeome.net/database/hccdb/download/patient/HCCDB18.patient.zip", pat_zip, mode = "wb", quiet = TRUE)

    unzip(expr_zip, exdir = raw_dir, overwrite = TRUE)
    unzip(samp_zip, exdir = raw_dir, overwrite = TRUE)
    unzip(pat_zip, exdir = raw_dir, overwrite = TRUE)
  }

  samp <- fread(samp_txt, data.table = FALSE)
  rownames(samp) <- samp[[1]]
  samp[[1]] <- NULL
  if (!("TYPE" %in% rownames(samp)) || !("PATIENT_ID" %in% rownames(samp))) {
    stop("[HCCDB18] sample 文件缺少 TYPE 或 PATIENT_ID 行，无法映射肿瘤样本到患者")
  }
  sample_type <- trimws(as.character(unlist(samp["TYPE", , drop = FALSE])))
  names(sample_type) <- colnames(samp)
  patient_id_by_sample <- as.character(unlist(samp["PATIENT_ID", , drop = FALSE]))
  names(patient_id_by_sample) <- colnames(samp)

  tumor_samples <- names(sample_type)[sample_type == "HCC"]
  if (length(tumor_samples) < 10) stop("[HCCDB18] 肿瘤样本数过少，无法用于外部验证")

  pat <- fread(pat_txt, data.table = FALSE)
  rownames(pat) <- pat[[1]]
  pat[[1]] <- NULL
  if (!("SUR" %in% rownames(pat)) || !("STATUS" %in% rownames(pat))) {
    stop("[HCCDB18] patient 文件缺少 SUR 或 STATUS 行，无法构建生存结局")
  }
  surv_time <- suppressWarnings(as.numeric(unlist(pat["SUR", , drop = FALSE])))
  names(surv_time) <- colnames(pat)
  status_text <- as.character(unlist(pat["STATUS", , drop = FALSE]))
  names(status_text) <- colnames(pat)
  clin <- data.frame(
    patient_id = names(surv_time),
    time_months = surv_time,
    status = ifelse(tolower(status_text) == "dead", 1, 0),
    stringsAsFactors = FALSE
  )
  rownames(clin) <- clin$patient_id

  expr_dt <- fread(expr_txt, data.table = FALSE)
  if (!all(c("Entrez_ID", "Symbol") %in% colnames(expr_dt))) {
    stop("[HCCDB18] 表达矩阵缺少 Entrez_ID/Symbol 列，无法解析")
  }
  genes <- make.unique(as.character(expr_dt$Symbol))
  expr_mat <- as.matrix(expr_dt[, -(1:2), drop = FALSE])
  rownames(expr_mat) <- genes

  tumor_samples <- intersect(tumor_samples, colnames(expr_mat))
  if (length(tumor_samples) < 10) stop("[HCCDB18] 肿瘤样本与表达矩阵交集过少，无法用于外部验证")

  expr_mat <- expr_mat[, tumor_samples, drop = FALSE]
  pid <- as.character(patient_id_by_sample[tumor_samples])
  keep <- !is.na(pid) & pid != ""
  expr_mat <- expr_mat[, keep, drop = FALSE]
  pid <- pid[keep]

  expr_t <- t(expr_mat)
  expr_sum <- rowsum(expr_t, pid)
  expr_avg <- expr_sum / as.vector(table(pid))
  expr_avg <- t(expr_avg)

  saveRDS(expr_avg, file.path(proc_dir, "HCCDB18_expr_symbol_hcc.rds"))
  saveRDS(clin, file.path(proc_dir, "HCCDB18_clinical_hcc.rds"))
}

filter_hcc_samples <- function(pdat) {
  text <- tolower(collapse_characteristics(pdat))
  is_hcc <- grepl("hcc|hepatocellular|hepatocarcinoma|tumou?r|cancer", text)
  exclude <- grepl("non[- ]?tumor|nontumor|normal|hepatitis|cirrhosis|fibrosis", text)
  keep <- is_hcc & !exclude
  keep
}

summarize_validation <- function(df, label) {
  df <- df %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))

  df$risk_group <- ifelse(df$risk_score >= median(df$risk_score), "High", "Low")

  surv_obj <- Surv(df$time_months, df$status)
  cox_fit <- coxph(surv_obj ~ risk_score, data = df)
  cindex <- summary(cox_fit)$concordance[1]

  logrank <- survdiff(surv_obj ~ risk_group, data = df)
  logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)

  max_time <- max(df$time_months, na.rm = TRUE)
  roc_times <- c(12, 36, 60)
  roc_times <- roc_times[roc_times < max_time]
  # Guard against unstable AUC estimates when very few samples remain at risk
  # (or too few events occur) at a given timepoint.
  min_at_risk <- 20
  min_events_by_t <- 10
  roc_times <- roc_times[vapply(roc_times, function(t) {
    at_risk <- sum(df$time_months >= t, na.rm = TRUE)
    events_by_t <- sum(df$status == 1 & df$time_months <= t, na.rm = TRUE)
    (at_risk >= min_at_risk) && (events_by_t >= min_events_by_t)
  }, logical(1))]

  auc_1y <- NA
  auc_3y <- NA
  auc_5y <- NA
  roc_obj <- NULL

  if (length(roc_times) > 0) {
    roc_obj <- timeROC(
      T = df$time_months,
      delta = df$status,
      marker = df$risk_score,
      cause = 1,
      weighting = "marginal",
      times = roc_times,
      iid = FALSE
    )

    if (12 %in% roc_times) auc_1y <- roc_obj$AUC[roc_obj$times == 12]
    if (36 %in% roc_times) auc_3y <- roc_obj$AUC[roc_obj$times == 36]
    if (60 %in% roc_times) auc_5y <- roc_obj$AUC[roc_obj$times == 60]
  }

  stats <- data.frame(
    Dataset = label,
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    Genes_Used_N = if ("Genes_Used_N" %in% colnames(df)) df$Genes_Used_N[1] else NA_integer_,
    Missing_Model_Genes = if ("Missing_Model_Genes" %in% colnames(df)) df$Missing_Model_Genes[1] else NA_character_,
    AUC_1y = auc_1y,
    AUC_3y = auc_3y,
    AUC_5y = auc_5y,
    Logrank_P = logrank_p,
    C_index = cindex
  )

  list(stats = stats, data = df, roc = roc_obj)
}

summarize_validation_dualcut <- function(df, label, train_cutoff) {
  df <- df %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))

  df$risk_group_median <- ifelse(df$risk_score >= median(df$risk_score), "High", "Low")
  df$risk_group_train <- ifelse(df$risk_score >= train_cutoff, "High", "Low")

  surv_obj <- Surv(df$time_months, df$status)
  cox_fit <- coxph(surv_obj ~ risk_score, data = df)
  cindex <- summary(cox_fit)$concordance[1]

  logrank_median <- survdiff(surv_obj ~ risk_group_median, data = df)
  logrank_p_median <- 1 - pchisq(logrank_median$chisq, length(logrank_median$n) - 1)

  logrank_p_train <- NA_real_
  if (length(unique(df$risk_group_train)) > 1) {
    logrank_train <- survdiff(surv_obj ~ risk_group_train, data = df)
    logrank_p_train <- 1 - pchisq(logrank_train$chisq, length(logrank_train$n) - 1)
  }

  max_time <- max(df$time_months, na.rm = TRUE)
  roc_times <- c(12, 36, 60)
  roc_times <- roc_times[roc_times < max_time]
  # Guard against unstable AUC estimates when very few samples remain at risk
  # (or too few events occur) at a given timepoint.
  min_at_risk <- 20
  min_events_by_t <- 10
  roc_times <- roc_times[vapply(roc_times, function(t) {
    at_risk <- sum(df$time_months >= t, na.rm = TRUE)
    events_by_t <- sum(df$status == 1 & df$time_months <= t, na.rm = TRUE)
    (at_risk >= min_at_risk) && (events_by_t >= min_events_by_t)
  }, logical(1))]

  auc_1y <- NA
  auc_3y <- NA
  auc_5y <- NA
  roc_obj <- NULL

  if (length(roc_times) > 0) {
    roc_obj <- timeROC(
      T = df$time_months,
      delta = df$status,
      marker = df$risk_score,
      cause = 1,
      weighting = "marginal",
      times = roc_times,
      iid = FALSE
    )

    if (12 %in% roc_times) auc_1y <- roc_obj$AUC[roc_obj$times == 12]
    if (36 %in% roc_times) auc_3y <- roc_obj$AUC[roc_obj$times == 36]
    if (60 %in% roc_times) auc_5y <- roc_obj$AUC[roc_obj$times == 60]
  }

  stats <- data.frame(
    Dataset = label,
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    Genes_Used_N = if ("Genes_Used_N" %in% colnames(df)) df$Genes_Used_N[1] else NA_integer_,
    Missing_Model_Genes = if ("Missing_Model_Genes" %in% colnames(df)) df$Missing_Model_Genes[1] else NA_character_,
    AUC_1y = auc_1y,
    AUC_3y = auc_3y,
    AUC_5y = auc_5y,
    Logrank_P_median = logrank_p_median,
    Logrank_P_traincut = logrank_p_train,
    C_index = cindex
  )

  list(stats = stats, data = df, roc = roc_obj)
}

plot_km <- function(df, title) {
  km_fit <- survfit(Surv(time_months, status) ~ risk_group, data = df)
  ggsurvplot(
    km_fit,
    data = df,
    pval = TRUE,
    risk.table = FALSE,
    palette = c("#2E9FDF", "#E7B800"),
    xlab = "Time (months)",
    ylab = "Overall Survival",
    legend.title = "",
    legend.labs = c("Low Risk", "High Risk"),
    ggtheme = theme_bw(base_size = 10)
  )$plot +
    labs(title = title)
}

plot_roc <- function(df, roc_obj, title) {
  if (is.null(roc_obj)) {
    return(ggplot() + theme_void() + labs(title = paste0(title, " (ROC NA)")))
  }
  roc_df <- data.frame()
  for (t in roc_obj$times) {
    idx <- which(roc_obj$times == t)
    tmp <- data.frame(
      FPR = roc_obj$FP[, idx],
      TPR = roc_obj$TP[, idx],
      Time = paste0(t / 12, "-year (AUC=", round(roc_obj$AUC[idx], 3), ")")
    )
    roc_df <- rbind(roc_df, tmp)
  }

  ggplot(roc_df, aes(x = FPR, y = TPR, color = Time)) +
    geom_path(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity", color = "") +
    theme_bw(base_size = 10) +
    theme(legend.position = c(0.65, 0.25))
}

message("[外部验证] 开始 TCGA-LIHC 验证...")
expr_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_expr.rds"))
clin_tcga <- readRDS(file.path(proc_dir, "TCGA_LIHC_clinical.rds"))

# Frozen-scaling sensitivity analysis (train-mean/train-sd; log2(x+1) if needed)
train_cutoff <- median(read.csv(file.path(res_dir, "risk_score_data.csv"))$risk_score, na.rm = TRUE)
scaler <- train_frozen_scaler(model_genes)

patient_id <- substr(colnames(expr_tcga), 1, 12)
expr_tcga_t <- t(expr_tcga)
expr_tcga_sum <- rowsum(expr_tcga_t, patient_id)
expr_tcga_avg <- expr_tcga_sum / as.vector(table(patient_id))
expr_tcga_avg <- t(expr_tcga_avg)

risk_tcga <- calc_risk_score(zscore_by_gene(expr_tcga_avg), coef_vec, cohort_label = "TCGA-LIHC")
tcga_df <- data.frame(
  patient_id = names(risk_tcga),
  risk_score = as.numeric(risk_tcga)
)
tcga_df <- add_gene_coverage(tcga_df, risk_tcga)

clin_tcga$patient_id <- clin_tcga$bcr_patient_barcode
clin_tcga$time_days <- ifelse(
  !is.na(clin_tcga$days_to_death),
  clin_tcga$days_to_death,
  clin_tcga$days_to_last_follow_up
)
clin_tcga$status <- ifelse(tolower(clin_tcga$vital_status) == "dead", 1, 0)
clin_tcga$time_months <- as.numeric(clin_tcga$time_days) / 30.4

val_tcga <- merge(tcga_df, clin_tcga, by = "patient_id")
val_tcga <- val_tcga %>% filter(!is.na(time_months), !is.na(status))

write.csv(val_tcga, file.path(res_dir, "TCGA_LIHC_risk_score.csv"), row.names = FALSE)

res_tcga <- summarize_validation(val_tcga, "TCGA-LIHC")

risk_tcga_frozen <- calc_risk_score(zscore_by_train(expr_tcga_avg, scaler), coef_vec, cohort_label = "TCGA-LIHC frozen")
tcga_df_frozen <- data.frame(
  patient_id = names(risk_tcga_frozen),
  risk_score = as.numeric(risk_tcga_frozen)
)
tcga_df_frozen <- add_gene_coverage(tcga_df_frozen, risk_tcga_frozen)
val_tcga_frozen <- merge(tcga_df_frozen, clin_tcga, by = "patient_id") %>%
  filter(!is.na(time_months), !is.na(status))
res_tcga_frozen <- summarize_validation_dualcut(val_tcga_frozen, "TCGA-LIHC", train_cutoff)

message("[外部验证] 开始 GSE76427 验证...")
expr_76427_path <- file.path(proc_dir, "GSE76427_expr_symbol.rds")
clin_76427_path <- file.path(proc_dir, "GSE76427_clinical.rds")

if (!file.exists(expr_76427_path) || !file.exists(clin_76427_path)) {
  g <- getGEO("GSE76427", GSEMatrix = TRUE, AnnotGPL = TRUE)
  es <- if (is.list(g)) g[[1]] else g
  pheno <- pData(es)
  expr_avg <- map_expression_to_symbols(es)

  pheno$time_months <- as.numeric(pheno[["duryears_os:ch1"]]) * 12
  pheno$status <- as.numeric(pheno[["event_os:ch1"]])

  saveRDS(expr_avg, expr_76427_path)
  saveRDS(pheno, clin_76427_path)
}

expr_76427 <- readRDS(expr_76427_path)
clin_76427 <- readRDS(clin_76427_path)

common_samples <- intersect(colnames(expr_76427), rownames(clin_76427))
expr_76427 <- expr_76427[, common_samples, drop = FALSE]
clin_76427 <- clin_76427[common_samples, , drop = FALSE]

risk_76427 <- calc_risk_score(zscore_by_gene(expr_76427), coef_vec, cohort_label = "GSE76427")
val_76427 <- data.frame(
  sample_id = names(risk_76427),
  risk_score = as.numeric(risk_76427),
  time_months = as.numeric(clin_76427$time_months),
  status = as.numeric(clin_76427$status)
)
val_76427 <- add_gene_coverage(val_76427, risk_76427)
val_76427 <- val_76427 %>% filter(!is.na(time_months), !is.na(status))

write.csv(val_76427, file.path(res_dir, "GSE76427_risk_score.csv"), row.names = FALSE)

res_76427 <- summarize_validation(val_76427, "GSE76427")

risk_76427_frozen <- calc_risk_score(zscore_by_train(expr_76427, scaler), coef_vec, cohort_label = "GSE76427 frozen")
val_76427_frozen <- data.frame(
  sample_id = names(risk_76427_frozen),
  risk_score = as.numeric(risk_76427_frozen),
  time_months = as.numeric(clin_76427$time_months[match(names(risk_76427_frozen), rownames(clin_76427))]),
  status = as.numeric(clin_76427$status[match(names(risk_76427_frozen), rownames(clin_76427))])
)
val_76427_frozen <- add_gene_coverage(val_76427_frozen, risk_76427_frozen)
val_76427_frozen <- val_76427_frozen %>% filter(!is.na(time_months), !is.na(status))
res_76427_frozen <- summarize_validation_dualcut(val_76427_frozen, "GSE76427", train_cutoff)

message("[外部验证] 开始 GSE10143 验证...")
expr_10143_path <- file.path(proc_dir, "GSE10143_expr_symbol.rds")
clin_10143_path <- file.path(proc_dir, "GSE10143_clinical.rds")
expr_10143_hcc_path <- file.path(proc_dir, "GSE10143_expr_symbol_hcc.rds")
clin_10143_hcc_path <- file.path(proc_dir, "GSE10143_clinical_hcc.rds")

expr_10143 <- readRDS(expr_10143_path)
clin_10143 <- readRDS(clin_10143_path)

keep_hcc <- filter_hcc_samples(clin_10143)
message("[外部验证] GSE10143 HCC 样本数: ", sum(keep_hcc), "/", nrow(clin_10143))
clin_10143 <- clin_10143[keep_hcc, , drop = FALSE]

common_samples <- intersect(colnames(expr_10143), rownames(clin_10143))
expr_10143 <- expr_10143[, common_samples, drop = FALSE]
clin_10143 <- clin_10143[common_samples, , drop = FALSE]

saveRDS(expr_10143, expr_10143_hcc_path)
saveRDS(clin_10143, clin_10143_hcc_path)

risk_10143 <- calc_risk_score(zscore_by_gene(expr_10143), coef_vec, cohort_label = "GSE10143-HCC")
val_10143 <- data.frame(
  sample_id = names(risk_10143),
  risk_score = as.numeric(risk_10143),
  time_months = as.numeric(clin_10143$time_months),
  status = as.numeric(clin_10143$status)
)
val_10143 <- add_gene_coverage(val_10143, risk_10143)
val_10143 <- val_10143 %>% filter(!is.na(time_months), !is.na(status))

write.csv(val_10143, file.path(res_dir, "GSE10143_HCC_risk_score.csv"), row.names = FALSE)

res_10143 <- summarize_validation(val_10143, "GSE10143-HCC")

risk_10143_frozen <- calc_risk_score(zscore_by_train(expr_10143, scaler), coef_vec, cohort_label = "GSE10143-HCC frozen")
val_10143_frozen <- data.frame(
  sample_id = names(risk_10143_frozen),
  risk_score = as.numeric(risk_10143_frozen),
  time_months = as.numeric(clin_10143$time_months[match(names(risk_10143_frozen), rownames(clin_10143))]),
  status = as.numeric(clin_10143$status[match(names(risk_10143_frozen), rownames(clin_10143))])
)
val_10143_frozen <- add_gene_coverage(val_10143_frozen, risk_10143_frozen)
val_10143_frozen <- val_10143_frozen %>% filter(!is.na(time_months), !is.na(status))
res_10143_frozen <- summarize_validation_dualcut(val_10143_frozen, "GSE10143-HCC", train_cutoff)

# Optional external cohort: GSE27150 (two-channel microarray; log2(Red/Green) ratio)
message("[外部验证] 开始 GSE27150 验证 (two-channel log2 ratio)...")
expr_27150_path <- file.path(proc_dir, "GSE27150_expr_symbol.rds")
clin_27150_path <- file.path(proc_dir, "GSE27150_clinical.rds")

if (!file.exists(expr_27150_path) || !file.exists(clin_27150_path)) {
  gse <- getGEO("GSE27150", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  expr_raw <- exprs(gse)
  # fData in GSE27150 matrix is empty; map via GPL
  expr_avg <- map_expression_with_gpl(expr_raw, "GPL13128")
  expr_avg <- maybe_log2p1(expr_avg)

  pheno <- pData(gse)
  # survival columns are already provided as "survival time:ch1" and "status (1-death,0-survival):ch1"
  if (!("survival time:ch1" %in% colnames(pheno)) || !("status (1-death,0-survival):ch1" %in% colnames(pheno))) {
    stop("GSE27150 缺少 survival time/status 字段，无法用于外部验证")
  }
  pheno$time_months <- suppressWarnings(as.numeric(pheno[["survival time:ch1"]]))
  pheno$status <- suppressWarnings(as.numeric(pheno[["status (1-death,0-survival):ch1"]]))
  rownames(pheno) <- rownames(pData(gse))

  saveRDS(expr_avg, expr_27150_path)
  saveRDS(pheno, clin_27150_path)
}

expr_27150 <- readRDS(expr_27150_path)
clin_27150 <- readRDS(clin_27150_path)

common_samples <- intersect(colnames(expr_27150), rownames(clin_27150))
expr_27150 <- expr_27150[, common_samples, drop = FALSE]
clin_27150 <- clin_27150[common_samples, , drop = FALSE]

risk_27150 <- calc_risk_score(zscore_by_gene(expr_27150), coef_vec, cohort_label = "GSE27150")
val_27150 <- data.frame(
  sample_id = names(risk_27150),
  risk_score = as.numeric(risk_27150),
  time_months = as.numeric(clin_27150$time_months[match(names(risk_27150), rownames(clin_27150))]),
  status = as.numeric(clin_27150$status[match(names(risk_27150), rownames(clin_27150))])
)
val_27150 <- add_gene_coverage(val_27150, risk_27150)
val_27150 <- val_27150 %>% filter(!is.na(time_months), !is.na(status))
write.csv(val_27150, file.path(res_dir, "GSE27150_risk_score.csv"), row.names = FALSE)
res_27150 <- summarize_validation(val_27150, "GSE27150")

risk_27150_frozen <- calc_risk_score(zscore_by_train(expr_27150, scaler), coef_vec, cohort_label = "GSE27150 frozen")
val_27150_frozen <- data.frame(
  sample_id = names(risk_27150_frozen),
  risk_score = as.numeric(risk_27150_frozen),
  time_months = as.numeric(clin_27150$time_months[match(names(risk_27150_frozen), rownames(clin_27150))]),
  status = as.numeric(clin_27150$status[match(names(risk_27150_frozen), rownames(clin_27150))])
)
val_27150_frozen <- add_gene_coverage(val_27150_frozen, risk_27150_frozen)
val_27150_frozen <- val_27150_frozen %>% filter(!is.na(time_months), !is.na(status))
res_27150_frozen <- summarize_validation_dualcut(val_27150_frozen, "GSE27150", train_cutoff)

message("[外部验证] 开始 ICGC-LIRI-JP (HCCDB18) RNA-seq 验证...")
expr_hccdb18_path <- file.path(proc_dir, "HCCDB18_expr_symbol_hcc.rds")
clin_hccdb18_path <- file.path(proc_dir, "HCCDB18_clinical_hcc.rds")
if (!file.exists(expr_hccdb18_path) || !file.exists(clin_hccdb18_path)) {
  prepare_hccdb18_lirijp()
}
expr_hccdb18 <- readRDS(expr_hccdb18_path)
clin_hccdb18 <- readRDS(clin_hccdb18_path)

common_patients <- intersect(colnames(expr_hccdb18), rownames(clin_hccdb18))
expr_hccdb18 <- expr_hccdb18[, common_patients, drop = FALSE]
clin_hccdb18 <- clin_hccdb18[common_patients, , drop = FALSE]

risk_hccdb18 <- calc_risk_score(zscore_by_gene(expr_hccdb18), coef_vec, cohort_label = "ICGC-LIRI-JP (HCCDB18)")
val_hccdb18 <- data.frame(
  patient_id = names(risk_hccdb18),
  risk_score = as.numeric(risk_hccdb18),
  time_months = as.numeric(clin_hccdb18$time_months[match(names(risk_hccdb18), rownames(clin_hccdb18))]),
  status = as.numeric(clin_hccdb18$status[match(names(risk_hccdb18), rownames(clin_hccdb18))])
)
val_hccdb18 <- add_gene_coverage(val_hccdb18, risk_hccdb18)
val_hccdb18 <- val_hccdb18 %>% filter(!is.na(time_months), !is.na(status))
write.csv(val_hccdb18, file.path(res_dir, "HCCDB18_LIRIJP_risk_score.csv"), row.names = FALSE)
res_hccdb18 <- summarize_validation(val_hccdb18, "ICGC-LIRI-JP (HCCDB18)")

risk_hccdb18_frozen <- calc_risk_score(zscore_by_train(expr_hccdb18, scaler), coef_vec, cohort_label = "ICGC-LIRI-JP (HCCDB18) frozen")
val_hccdb18_frozen <- data.frame(
  patient_id = names(risk_hccdb18_frozen),
  risk_score = as.numeric(risk_hccdb18_frozen),
  time_months = as.numeric(clin_hccdb18$time_months[match(names(risk_hccdb18_frozen), rownames(clin_hccdb18))]),
  status = as.numeric(clin_hccdb18$status[match(names(risk_hccdb18_frozen), rownames(clin_hccdb18))])
)
val_hccdb18_frozen <- add_gene_coverage(val_hccdb18_frozen, risk_hccdb18_frozen)
val_hccdb18_frozen <- val_hccdb18_frozen %>% filter(!is.na(time_months), !is.na(status))
res_hccdb18_frozen <- summarize_validation_dualcut(val_hccdb18_frozen, "ICGC-LIRI-JP (HCCDB18)", train_cutoff)

stats <- bind_rows(res_tcga$stats, res_76427$stats, res_10143$stats, res_27150$stats, res_hccdb18$stats)
write.csv(stats, file.path(res_dir, "external_validation_stats.csv"), row.names = FALSE)

stats_frozen <- bind_rows(
  transform(res_tcga_frozen$stats, Scaling = "Train-frozen"),
  transform(res_76427_frozen$stats, Scaling = "Train-frozen"),
  transform(res_10143_frozen$stats, Scaling = "Train-frozen"),
  transform(res_27150_frozen$stats, Scaling = "Train-frozen"),
  transform(res_hccdb18_frozen$stats, Scaling = "Train-frozen")
)
write.csv(stats_frozen, file.path(res_dir, "external_validation_stats_frozen_scaling.csv"), row.names = FALSE)
write.csv(stats_frozen, file.path(res_dir, "Supplementary_Table_ExternalValidation_FrozenScaling.csv"), row.names = FALSE)

p_tcga_km <- plot_km(res_tcga$data, "TCGA-LIHC")
p_tcga_roc <- plot_roc(res_tcga$data, res_tcga$roc, "TCGA-LIHC ROC")

p_76427_km <- plot_km(res_76427$data, "GSE76427")
p_76427_roc <- plot_roc(res_76427$data, res_76427$roc, "GSE76427 ROC")

p_10143_km <- plot_km(res_10143$data, "GSE10143-HCC")
p_10143_roc <- plot_roc(res_10143$data, res_10143$roc, "GSE10143-HCC ROC")

p_27150_km <- plot_km(res_27150$data, "GSE27150 (log2 ratio)")
p_27150_roc <- plot_roc(res_27150$data, res_27150$roc, "GSE27150 ROC")

p_hccdb18_km <- plot_km(res_hccdb18$data, "ICGC-LIRI-JP (HCCDB18)")
p_hccdb18_roc <- plot_roc(res_hccdb18$data, res_hccdb18$roc, "ICGC-LIRI-JP ROC")

fig <- plot_grid(
  p_tcga_km, p_tcga_roc,
  p_76427_km, p_76427_roc,
  p_10143_km, p_10143_roc,
  p_27150_km, p_27150_roc,
  p_hccdb18_km, p_hccdb18_roc,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
)

ggsave(file.path(plot_dir, "FigureS_external_validation.pdf"), fig, width = 12, height = 21)
ggsave(file.path(plot_dir, "FigureS_external_validation.png"), fig, width = 12, height = 21, dpi = 300, bg = "white")

message("[外部验证] 完成")
message("  ✅ results/TCGA_LIHC_risk_score.csv")
message("  ✅ results/GSE76427_risk_score.csv")
message("  ✅ results/GSE10143_HCC_risk_score.csv")
message("  ✅ results/GSE27150_risk_score.csv")
message("  ✅ results/HCCDB18_LIRIJP_risk_score.csv")
message("  ✅ results/external_validation_stats.csv")
message("  ✅ plots/supplementary/FigureS_external_validation.pdf")
message("  ✅ plots/supplementary/FigureS_external_validation.png")
