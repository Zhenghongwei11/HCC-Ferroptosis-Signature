#!/usr/bin/env Rscript

# 00c_prepare_hcc_ferro_genes.R
# Build HCC-specific ferroptosis gene list from FerrDb + GSE14520 DEG

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(AnnotationDbi)
  library(hgu133a2.db)
})

if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

raw_dir <- "data/raw"
proc_dir <- "data/processed"
ref_dir <- "data/references"
res_dir <- "results"

ferrdb_dir <- file.path(raw_dir, "ferrdb_v2_early_preview_20231231")
ferro_ref_file <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")

# 1) Read FerrDb-derived gene list. Public releases use the curated
# data/references/ferroptosis_genes_expanded.csv table to avoid bundling
# large raw database annotation text that is not needed for reproduction.
read_symbols <- function(path, symbol_col) {
  if (!file.exists(path)) return(character())
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (!symbol_col %in% colnames(df)) return(character())
  sym <- df[[symbol_col]]
  sym <- trimws(sym)
  sym <- sym[!is.na(sym) & sym != ""]
  toupper(sym)
}

if (file.exists(ferro_ref_file)) {
  ferrdb_ref <- read.csv(ferro_ref_file, stringsAsFactors = FALSE)
  if (!("Gene" %in% colnames(ferrdb_ref))) stop("ferroptosis_genes_expanded.csv 缺少 Gene 列")
  ferrdb_genes <- unique(toupper(trimws(ferrdb_ref$Gene)))
  ferrdb_genes <- ferrdb_genes[!is.na(ferrdb_genes) & ferrdb_genes != ""]
} else if (dir.exists(ferrdb_dir)) {
  ferrdb_genes <- unique(c(
    read_symbols(file.path(ferrdb_dir, "driver.csv"), "Symbol_or_reported_abbr"),
    read_symbols(file.path(ferrdb_dir, "suppressor.csv"), "Symbol"),
    read_symbols(file.path(ferrdb_dir, "marker.csv"), "Symbol"),
    read_symbols(file.path(ferrdb_dir, "unclassified.reg.csv"), "Symbol")
  ))
} else {
  stop("缺少 ferroptosis reference: data/references/ferroptosis_genes_expanded.csv")
}

if (length(ferrdb_genes) == 0) {
  stop("FerrDb-derived gene list is empty")
}

# 2) Ensure DEG table exists; if not, compute from GSE14520
ensure_deg <- function() {
  deg_path <- file.path(res_dir, "deg_GSE14520_all.csv")
  if (file.exists(deg_path)) return(deg_path)

  expr_path <- file.path(proc_dir, "GSE14520_expr.rds")
  clin_path <- file.path(proc_dir, "GSE14520_clinical.rds")
  if (!file.exists(expr_path) || !file.exists(clin_path)) {
    stop("缺少 GSE14520 处理后文件，无法计算 DEG")
  }

  expr <- readRDS(expr_path)
  pdata <- readRDS(clin_path)

  group <- ifelse(
    grepl("Non-Tumor|non-tumor|Normal|normal|control|healthy", pdata$title, ignore.case = TRUE),
    "Normal",
    "HCC"
  )
  group <- factor(group, levels = c("Normal", "HCC"))

  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(expr, design)
  contrast <- makeContrasts(HCC - Normal, levels = design)
  fit <- eBayes(contrasts.fit(fit, contrast))

  deg <- topTable(fit, adjust.method = "fdr", number = Inf)
  deg$ProbeID <- rownames(deg)
  deg$Gene <- mapIds(
    hgu133a2.db,
    keys = rownames(deg),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  deg <- deg %>% filter(!is.na(Gene) & Gene != "")
  write.csv(deg, deg_path, row.names = FALSE)
  deg_path
}

deg_path <- ensure_deg()

# 3) Filter ferroptosis genes by HCC DEG (adj.P < 0.05)
deg <- read.csv(deg_path, stringsAsFactors = FALSE)
deg$Gene <- toupper(deg$Gene)

deg_sig <- deg %>% dplyr::filter(adj.P.Val < 0.05)

hcc_ferro <- deg_sig %>%
  dplyr::filter(Gene %in% ferrdb_genes) %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_min(order_by = adj.P.Val, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    HCC_logFC = logFC,
    HCC_adjP = adj.P.Val,
    HCC_direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  dplyr::select(Gene, HCC_logFC, HCC_adjP, HCC_direction)

if (nrow(hcc_ferro) == 0) {
  stop("HCC 特异 ferroptosis 基因为空，请检查阈值或数据")
}

out_path <- file.path(res_dir, "ferroptosis_genes_hcc_context.csv")
write.csv(hcc_ferro, out_path, row.names = FALSE)

message("[00c] ✅ ferroptosis_genes_hcc_context.csv 已生成")
message("  记录数: ", nrow(hcc_ferro))
message("  输出: ", out_path)
