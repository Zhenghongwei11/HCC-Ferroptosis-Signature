#!/usr/bin/env Rscript

# 02b_multi_cohort_DEG.R
# 单队列 DEG 分析（GSE14520） + 与铁死亡基因交集
# 输入: GSE14520_expr.rds, GSE14520_clinical.rds
# 输出: deg_GSE14520_all.csv, deg_GSE14520_sig.csv, ferroptosis_DEG_intersection.csv

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(AnnotationDbi)
  library(hgu133a2.db)  # Affymetrix HG-U133A 2.0 注释包 (GPL571)
})

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir <- "results"
ref_dir <- "data/references"

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

message("[DEG] 开始加载数据...")

# 探针ID到基因符号的转换函数
probe_to_symbol <- function(probe_ids) {
  tryCatch({
    symbols <- mapIds(
      hgu133a2.db,
      keys = probe_ids,
      column = "SYMBOL",
      keytype = "PROBEID",
      multiVals = "first"
    )
    return(symbols)
  }, error = function(e) {
    message("  注释包转换失败，返回NA")
    return(rep(NA, length(probe_ids)))
  })
}

# 1. 加载表达矩阵和临床信息
expr_14 <- readRDS(file.path(proc_dir, "GSE14520_expr.rds"))
pdata_14 <- readRDS(file.path(proc_dir, "GSE14520_clinical.rds"))

# 2. 加载铁死亡基因集（优先使用 HCC-context；缺省回退到参考 FerrDb 列表）
ferro_context_path <- file.path(res_dir, "ferroptosis_genes_hcc_context.csv")
ferro_ref_path <- file.path(ref_dir, "ferroptosis_genes_expanded.csv")
ferro_path <- if (file.exists(ferro_context_path)) ferro_context_path else ferro_ref_path
if (!file.exists(ferro_path)) {
  stop("缺少 ferroptosis 基因集：请先运行 00c，或准备 data/references/ferroptosis_genes_expanded.csv")
}
ferro_genes <- read.csv(ferro_path)$Gene

message("[DEG] GSE14520: ", dim(expr_14)[1], " x ", dim(expr_14)[2])
message("[DEG] 铁死亡基因(去重): ", length(unique(ferro_genes)))

# ============================================
# GSE14520 DEG分析
# ============================================
message("[DEG] 开始GSE14520 DEG分析...")

# 构造分组信息（HCC vs Normal）
# GSE14520样本标记为 "Liver Tumor Tissue" 和 "Liver Non-Tumor Tissue"
group_14 <- ifelse(
  grepl("Non-Tumor|non-tumor|Normal|normal|control|healthy", pdata_14$title, ignore.case = TRUE),
  "Normal",
  "HCC"
)
group_14 <- factor(group_14, levels = c("Normal", "HCC"))

message("[DEG] GSE14520分组: Normal=", sum(group_14 == "Normal"), ", HCC=", sum(group_14 == "HCC"))

# Limma DEG分析
design_14 <- model.matrix(~0 + group_14)
colnames(design_14) <- levels(group_14)
fit_14 <- lmFit(expr_14, design_14)
contrast_14 <- makeContrasts(HCC - Normal, levels = design_14)
fit_14 <- eBayes(contrasts.fit(fit_14, contrast_14))

deg_14 <- topTable(fit_14, adjust.method = "fdr", number = Inf)

# 添加基因符号列
message("[DEG] 转换GSE14520探针ID到基因符号...")
deg_14$ProbeID <- rownames(deg_14)
deg_14$Gene <- probe_to_symbol(rownames(deg_14))
deg_14 <- deg_14 %>% filter(!is.na(Gene) & Gene != "")

deg_14_sig <- deg_14 %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)

message("[DEG] GSE14520显著DEG: ", nrow(deg_14_sig))

write.csv(deg_14, file.path(res_dir, "deg_GSE14520_all.csv"), row.names = FALSE)
write.csv(deg_14_sig, file.path(res_dir, "deg_GSE14520_sig.csv"), row.names = FALSE)

# ============================================
# 与铁死亡基因交集
# ============================================
message("[DEG] 进行与铁死亡基因的交集...")

ferro_deg_14 <- intersect(deg_14_sig$Gene, ferro_genes)
message("[DEG] GSE14520铁死亡DEG: ", length(ferro_deg_14))

ferro_deg_result <- data.frame(
  Gene = ferro_deg_14,
  logFC_14 = deg_14_sig$logFC[match(ferro_deg_14, deg_14_sig$Gene)],
  adj.P.Val_14 = deg_14_sig$adj.P.Val[match(ferro_deg_14, deg_14_sig$Gene)]
)

write.csv(ferro_deg_result, file.path(res_dir, "ferroptosis_DEG_intersection.csv"), row.names = FALSE)

message("[DEG] 完成")
message("  ✅ ", file.path(res_dir, "deg_GSE14520_all.csv"))
message("  ✅ ", file.path(res_dir, "deg_GSE14520_sig.csv"))
message("  ✅ ", file.path(res_dir, "ferroptosis_DEG_intersection.csv"))
