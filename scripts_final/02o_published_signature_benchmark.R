#!/usr/bin/env Rscript

# 02o_published_signature_benchmark.R
# Head-to-head benchmark against published HCC prognostic signatures.
# Published formulas are evaluated with the same cohort-level gene
# standardization and survival endpoints used for the current validation set.

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(timeROC)
  library(ggplot2)
})

proc_dir <- "data/processed"
res_dir <- "results"
plot_dir <- "plots/supplementary"
bench_dir <- file.path(res_dir, "benchmarks")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bench_dir, showWarnings = FALSE, recursive = TRUE)

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

zscore_by_gene <- function(expr_mat) {
  z <- t(scale(t(expr_mat)))
  z[is.na(z)] <- 0
  z
}

standardize_gene_symbols <- function(expr_mat) {
  rownames(expr_mat) <- toupper(rownames(expr_mat))
  expr_mat
}

collapse_duplicate_genes <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  if (anyDuplicated(rownames(expr_mat)) == 0) return(expr_mat)
  expr_sum <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)
  expr_cnt <- rowsum(matrix(1, nrow(expr_mat), ncol(expr_mat)), group = rownames(expr_mat), reorder = FALSE)
  expr_sum / expr_cnt
}

apply_gene_aliases <- function(expr_mat) {
  aliases <- list(
    CARS1 = c("CARS1", "CARS"),
    HILPDA = c("HILPDA", "HIG2", "C7ORF68"),
    PKM = c("PKM", "PKM2")
  )
  for (target in names(aliases)) {
    if (target %in% rownames(expr_mat)) next
    hit <- aliases[[target]][aliases[[target]] %in% rownames(expr_mat)]
    hit <- setdiff(hit, target)
    if (length(hit) > 0) {
      expr_mat <- rbind(expr_mat, expr_mat[hit[1], , drop = FALSE])
      rownames(expr_mat)[nrow(expr_mat)] <- target
    }
  }
  collapse_duplicate_genes(expr_mat)
}

prep_expr <- function(expr_mat) {
  expr_mat %>%
    standardize_gene_symbols() %>%
    collapse_duplicate_genes() %>%
    apply_gene_aliases()
}

read_risk_csv <- function(path, dataset) {
  df <- safe_read_csv(path)
  if ("time" %in% names(df) && !"time_months" %in% names(df)) df$time_months <- df$time
  id_col <- intersect(c("sample", "sample_id", "patient_id"), names(df))[1]
  if (is.na(id_col)) stop("Risk file lacks a sample or patient identifier: ", path)
  if (!all(c("risk_score", "time_months", "status") %in% names(df))) {
    stop("Risk file lacks risk_score/time_months/status columns: ", path)
  }
  df %>%
    transmute(
      Dataset = dataset,
      sample_id = .data[[id_col]],
      time_months = as.numeric(time_months),
      status = as.integer(status),
      score = as.numeric(risk_score)
    ) %>%
    filter(!is.na(time_months), !is.na(status), !is.na(score))
}

tcga_patient_expr <- function(expr_mat) {
  patient_id <- substr(colnames(expr_mat), 1, 12)
  expr_t <- t(expr_mat)
  expr_sum <- rowsum(expr_t, patient_id)
  t(expr_sum / as.vector(table(patient_id)))
}

load_cohort_expr <- function(dataset) {
  if (dataset == "GSE14520") {
    prep_expr(readRDS(file.path(proc_dir, "GSE14520_expr_symbol.rds")))
  } else if (dataset == "TCGA-LIHC") {
    prep_expr(tcga_patient_expr(readRDS(file.path(proc_dir, "TCGA_LIHC_expr.rds"))))
  } else if (dataset == "GSE76427") {
    prep_expr(readRDS(file.path(proc_dir, "GSE76427_expr_symbol.rds")))
  } else if (dataset == "GSE10143-HCC") {
    prep_expr(readRDS(file.path(proc_dir, "GSE10143_expr_symbol_hcc.rds")))
  } else if (dataset == "GSE27150") {
    prep_expr(readRDS(file.path(proc_dir, "GSE27150_expr_symbol.rds")))
  } else if (dataset == "ICGC-LIRI-JP (HCCDB18)") {
    prep_expr(readRDS(file.path(proc_dir, "HCCDB18_expr_symbol_hcc.rds")))
  } else {
    stop("Unknown dataset: ", dataset)
  }
}

published_signatures <- list(
  "Current nine-gene score" = list(
    class = "Current model",
    citation = "This study",
    doi = NA_character_,
    coefficients = NA
  ),
  "Liang 2020 10-FRG" = list(
    class = "Published ferroptosis signature",
    citation = "Liang et al. 2020, Int J Biol Sci",
    doi = "10.7150/ijbs.45050",
    coefficients = c(
      SLC7A11 = 0.105, G6PD = 0.116, CISD1 = 0.106, CARS1 = 0.076,
      SLC1A5 = 0.077, ACACA = 0.092, ACSL3 = 0.005, NQO1 = 0.006,
      NFS1 = 0.087, GPX4 = 0.135
    )
  ),
  "Zhang 2022 5-FRG" = list(
    class = "Published ferroptosis signature",
    citation = "Zhang et al. 2022, Front Mol Biosci",
    doi = "10.3389/fmolb.2022.940575",
    coefficients = c(
      SLC7A11 = 0.1668, SLC1A5 = 0.1507, TFRC = 0.0221,
      CARS1 = 0.1515, RPL8 = 0.0234
    )
  ),
  "Xu 2022 13-FRG/PRG" = list(
    class = "Published ferroptosis/pyroptosis signature",
    citation = "Xu et al. 2022, BMC Cancer",
    doi = "10.1186/s12885-022-09301-0",
    coefficients = c(
      ATG3 = 0.0599818988151381,
      FLT3 = -0.321132320389413,
      G6PD = 0.0881814324303116,
      GLMN = 0.130781902193193,
      HILPDA = 0.119282064768739,
      LRPPRC = 0.00792886569542188,
      MKI67 = 0.0165502840606549,
      NRAS = 0.0916391974243284,
      PRDX6 = 0.114398632925529,
      SLC1A5 = 0.0521211560497305,
      SLC7A11 = 0.0616722848553423,
      SQSTM1 = 0.00940765304518399,
      UBE2D2 = 0.04686256927426
    )
  ),
  "Liu 2019 6-gene HCC" = list(
    class = "Published general HCC signature",
    citation = "Liu et al. 2019, Cancer Cell Int",
    doi = "10.1186/s12935-019-0858-2",
    coefficients = c(
      CSE1L = 0.0606, CSTB = 0.0257, MTHFR = 0.1177,
      DAGLA = 0.1912, MMP10 = 0.4324, GYS2 = -0.1003
    )
  ),
  "MKI67 single-gene baseline" = list(
    class = "Proliferation baseline",
    citation = "Single-gene proliferation marker",
    doi = NA_character_,
    coefficients = c(MKI67 = 1)
  ),
  "Cell-cycle metagene baseline" = list(
    class = "Proliferation baseline",
    citation = "Cell-cycle metagene baseline",
    doi = NA_character_,
    coefficients = c(
      MKI67 = 1, TOP2A = 1, PCNA = 1, CCNB1 = 1, CCNB2 = 1, CDK1 = 1,
      CDC20 = 1, BIRC5 = 1, AURKA = 1, AURKB = 1, CCNA2 = 1, MCM2 = 1, MCM6 = 1
    )
  )
)

cohort_files <- c(
  "GSE14520" = file.path(res_dir, "risk_score_data.csv"),
  "TCGA-LIHC" = file.path(res_dir, "TCGA_LIHC_risk_score.csv"),
  "GSE76427" = file.path(res_dir, "GSE76427_risk_score.csv"),
  "GSE10143-HCC" = file.path(res_dir, "GSE10143_HCC_risk_score.csv"),
  "GSE27150" = file.path(res_dir, "GSE27150_risk_score.csv"),
  "ICGC-LIRI-JP (HCCDB18)" = file.path(res_dir, "HCCDB18_LIRIJP_risk_score.csv")
)

score_signature <- function(expr_z, coef_vec) {
  coef_vec <- coef_vec[!is.na(names(coef_vec))]
  names(coef_vec) <- toupper(names(coef_vec))
  available <- intersect(names(coef_vec), rownames(expr_z))
  missing <- setdiff(names(coef_vec), rownames(expr_z))
  if (length(available) == 0) {
    return(list(score = NULL, available = available, missing = missing))
  }
  score <- colSums(expr_z[available, , drop = FALSE] * coef_vec[available])
  list(score = score, available = available, missing = missing)
}

safe_auc <- function(time_months, status, marker, t_months) {
  ok <- !is.na(time_months) & !is.na(status) & !is.na(marker)
  if (sum(ok) < 20) return(NA_real_)
  time_months <- time_months[ok]
  status <- status[ok]
  marker <- marker[ok]
  if (sum(time_months <= t_months & status == 1) < 5 || sum(time_months > t_months) < 5) return(NA_real_)
  out <- try(timeROC(
    T = time_months,
    delta = status,
    marker = marker,
    cause = 1,
    times = t_months,
    iid = FALSE
  ), silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  as.numeric(out$AUC[2])
}

score_metrics <- function(df) {
  df <- df %>% filter(!is.na(score), !is.na(time_months), !is.na(status))
  out <- tibble(
    N = nrow(df),
    Events = sum(df$status == 1),
    C_index = NA_real_,
    C_index_SE = NA_real_,
    HR_per_1SD = NA_real_,
    CI95_L = NA_real_,
    CI95_U = NA_real_,
    P_value = NA_real_,
    LogRank_P_median = NA_real_,
    AUC_1y = NA_real_,
    AUC_3y = NA_real_,
    AUC_5y = NA_real_
  )
  if (nrow(df) < 20 || length(unique(df$status)) < 2 || sd(df$score, na.rm = TRUE) == 0) return(out)
  df$score_z <- as.numeric(scale(df$score))
  fit <- try(coxph(Surv(time_months, status) ~ score_z, data = df), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    s <- summary(fit)
    out$C_index <- as.numeric(s$concordance[1])
    out$C_index_SE <- as.numeric(s$concordance[2])
    out$HR_per_1SD <- s$conf.int[1, "exp(coef)"]
    out$CI95_L <- s$conf.int[1, "lower .95"]
    out$CI95_U <- s$conf.int[1, "upper .95"]
    out$P_value <- s$coefficients[1, "Pr(>|z|)"]
  }
  if (nrow(df) >= 20 && length(unique(df$status)) == 2) {
    df$risk_group <- ifelse(df$score >= median(df$score, na.rm = TRUE), "High", "Low")
    lr <- try(survdiff(Surv(time_months, status) ~ risk_group, data = df), silent = TRUE)
    if (!inherits(lr, "try-error")) out$LogRank_P_median <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  }
  out$AUC_1y <- safe_auc(df$time_months, df$status, df$score, 12)
  out$AUC_3y <- safe_auc(df$time_months, df$status, df$score, 36)
  out$AUC_5y <- safe_auc(df$time_months, df$status, df$score, 60)
  out
}

source_rows <- bind_rows(lapply(names(published_signatures), function(model_name) {
  sig <- published_signatures[[model_name]]
  tibble(
    Model = model_name,
    Model_class = sig$class,
    Citation = sig$citation,
    DOI = sig$doi,
    Genes = if (all(is.na(sig$coefficients))) {
      paste(safe_read_csv(file.path(res_dir, "prognostic_model_coef.csv"))$Gene, collapse = ";")
    } else {
      paste(names(sig$coefficients), collapse = ";")
    },
    Coefficients = if (all(is.na(sig$coefficients))) {
      paste(
        safe_read_csv(file.path(res_dir, "prognostic_model_coef.csv"))$Coefficient,
        collapse = ";"
      )
    } else {
      paste(unname(sig$coefficients), collapse = ";")
    },
    Formula_source = case_when(
      Model == "Liang 2020 10-FRG" ~ "PMC7378635 text; linear predictor used because exponentiation is monotonic",
      Model == "Zhang 2022 5-FRG" ~ "Frontiers full text, TCGA model formula",
      Model == "Xu 2022 13-FRG/PRG" ~ "PMC8892773 text",
      Model == "Liu 2019 6-gene HCC" ~ "Springer full text",
      TRUE ~ "Current analysis definition"
    )
  )
}))
write.csv(source_rows, file.path(res_dir, "published_signature_sources.csv"), row.names = FALSE)

benchmark_rows <- list()
coverage_rows <- list()

for (dataset in names(cohort_files)) {
  message("[Signature benchmark] ", dataset)
  surv_df <- read_risk_csv(cohort_files[[dataset]], dataset)
  expr <- load_cohort_expr(dataset)
  common_samples <- intersect(surv_df$sample_id, colnames(expr))
  if (length(common_samples) < 20) stop("Too few expression/survival matches for ", dataset)
  surv_df <- surv_df %>% filter(sample_id %in% common_samples)
  expr <- expr[, surv_df$sample_id, drop = FALSE]
  expr_z <- zscore_by_gene(expr)

  for (model_name in names(published_signatures)) {
    sig <- published_signatures[[model_name]]
    if (model_name == "Current nine-gene score") {
      score_df <- surv_df
      total_genes <- strsplit(source_rows$Genes[source_rows$Model == model_name], ";", fixed = TRUE)[[1]]
      risk_source <- safe_read_csv(cohort_files[[dataset]])
      if ("Genes_Used" %in% names(risk_source)) {
        available <- strsplit(as.character(unique(risk_source$Genes_Used)[1]), ";", fixed = TRUE)[[1]]
      } else {
        available <- intersect(total_genes, rownames(expr_z))
      }
      if (length(available) == 1 && (is.na(available) || available == "")) available <- intersect(total_genes, rownames(expr_z))
      missing <- setdiff(total_genes, available)
    } else {
      scored <- score_signature(expr_z, sig$coefficients)
      available <- scored$available
      missing <- scored$missing
      total_genes <- names(sig$coefficients)
      if (is.null(scored$score)) next
      score_df <- surv_df
      score_df$score <- as.numeric(scored$score[score_df$sample_id])
    }

    metrics <- score_metrics(score_df)
    coverage <- length(available) / length(total_genes)
    benchmark_rows[[length(benchmark_rows) + 1]] <- bind_cols(
      tibble(
        Dataset = dataset,
        Model = model_name,
        Model_class = sig$class,
        Genes_Used_N = length(available),
        Genes_Total_N = length(total_genes),
        Gene_Coverage = coverage,
        Genes_Used = paste(available, collapse = ";"),
        Missing_Genes = paste(missing, collapse = ";")
      ),
      metrics
    )
    coverage_rows[[length(coverage_rows) + 1]] <- tibble(
      Dataset = dataset,
      Model = model_name,
      Model_class = sig$class,
      Genes_Used_N = length(available),
      Genes_Total_N = length(total_genes),
      Gene_Coverage = coverage,
      Genes_Used = paste(available, collapse = ";"),
      Missing_Genes = paste(missing, collapse = ";")
    )
  }
}

benchmark_df <- bind_rows(benchmark_rows) %>%
  mutate(
    Benchmark_Interpretation = case_when(
      Gene_Coverage >= 1 ~ "Complete formula coverage",
      Gene_Coverage >= 0.8 ~ "Near-complete formula coverage",
      Gene_Coverage >= 0.5 ~ "Partial formula coverage; interpret cautiously",
      TRUE ~ "Low formula coverage; sensitivity only"
    )
  )

coverage_df <- bind_rows(coverage_rows) %>%
  mutate(
    Benchmark_Interpretation = case_when(
      Gene_Coverage >= 1 ~ "Complete formula coverage",
      Gene_Coverage >= 0.8 ~ "Near-complete formula coverage",
      Gene_Coverage >= 0.5 ~ "Partial formula coverage; interpret cautiously",
      TRUE ~ "Low formula coverage; sensitivity only"
    )
  )

write.csv(benchmark_df, file.path(res_dir, "published_signature_benchmark.csv"), row.names = FALSE)
write.csv(coverage_df, file.path(res_dir, "published_signature_gene_coverage.csv"), row.names = FALSE)
write.table(
  benchmark_df,
  file.path(bench_dir, "published_signature_head_to_head.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

plot_df <- benchmark_df %>%
  filter(Gene_Coverage >= 0.8) %>%
  mutate(
    Dataset = factor(Dataset, levels = names(cohort_files)),
    Model = factor(Model, levels = names(published_signatures))
  )

p <- ggplot(plot_df, aes(x = Dataset, y = C_index, color = Model, shape = Model_class)) +
  geom_hline(yintercept = 0.5, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_point(size = 2.4, position = position_dodge(width = 0.55), na.rm = TRUE) +
  coord_cartesian(ylim = c(0.40, 0.78)) +
  scale_color_manual(values = c(
    "Current nine-gene score" = "#1F5A68",
    "Liang 2020 10-FRG" = "#7A5C9E",
    "Zhang 2022 5-FRG" = "#B0663B",
    "Xu 2022 13-FRG/PRG" = "#577A35",
    "Liu 2019 6-gene HCC" = "#A33F4A",
    "MKI67 single-gene baseline" = "#5F6C72",
    "Cell-cycle metagene baseline" = "#2F7D6D"
  )) +
  labs(x = NULL, y = "C-index", color = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "Supp_Figure_S13_Published_Signature_Benchmark.pdf"), p, width = 10.5, height = 6.5, device = cairo_pdf)
ggsave(file.path(plot_dir, "Supp_Figure_S13_Published_Signature_Benchmark.png"), p, width = 10.5, height = 6.5, dpi = 300, bg = "white")

message("[Signature benchmark] Published-signature benchmark finished.")
