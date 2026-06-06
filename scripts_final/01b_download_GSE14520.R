#!/usr/bin/env Rscript

# 01b_download_GSE14520.R
# 下载并整理 GSE14520（HBV-related HCC，大队列，含临床）的表达矩阵与临床数据
# 输出：
#   - data/raw/GSE14520/GSE14520_series_matrix.rds (可选)
#   - data/processed/GSE14520_expr.rds
#   - data/processed/GSE14520_clinical.rds
#   - data/processed/GSE14520_expr_log2.rds（如需）
#   - 一个简单的临床字段汇总打印，方便人工确认

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(stringr)
})

message("[GSE14520] 开始下载与解析...")

out_raw_dir <- "data/raw/GSE14520"
out_proc_dir <- "data/processed"
if (!dir.exists(out_raw_dir)) dir.create(out_raw_dir, recursive = TRUE)
if (!dir.exists(out_proc_dir)) dir.create(out_proc_dir, recursive = TRUE)

# 本地文件已存在，优先解析 GPL3921（主平台）
file_gpl3921 <- file.path(out_raw_dir, "GSE14520-GPL3921_series_matrix.txt.gz")
file_gpl571 <- file.path(out_raw_dir, "GSE14520-GPL571_series_matrix.txt.gz")
eset_list <- list()
if (file.exists(file_gpl3921)) {
  message("[GSE14520] 解析 GPL3921 (主平台)...")
  eset3921 <- getGEO(filename = file_gpl3921, GSEMatrix = TRUE, getGPL = FALSE)
  eset_list[["GPL3921"]] <- eset3921
}
if (file.exists(file_gpl571)) {
  message("[GSE14520] 解析 GPL571...")
  eset571 <- getGEO(filename = file_gpl571, GSEMatrix = TRUE, getGPL = FALSE)
  eset_list[["GPL571"]] <- eset571
}
if (length(eset_list) == 0) {
  stop("[GSE14520] 未找到本地 series matrix 文件。")
}
message("[GSE14520] 共解析 ", length(eset_list), " 个平台。")


# 收集表达矩阵与临床
expr_list <- list()
clin_list <- list()

for (i in seq_along(eset_list)) {
  eset <- eset_list[[i]]
  platform <- annotation(eset)
  cat(sprintf("[GSE14520] 发现平台: %s, 样本数: %d, 基因数: %d\n",
              platform, ncol(exprs(eset)), nrow(exprs(eset))))

  expr <- exprs(eset)
  pdata <- pData(eset)

  # 清洗列名：若为探针，尝试使用特征名
  fdata <- fData(eset)
  gene_col_candidates <- c("Gene.symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL")
  gene_col <- intersect(gene_col_candidates, colnames(fdata))
  if (length(gene_col) > 0) {
    gene_symbols <- as.character(fdata[[gene_col[1]]])
    # 合并重复基因：取探针平均值
    if (!all(is.na(gene_symbols))) {
      valid <- !is.na(gene_symbols) & gene_symbols != ""
      expr_valid <- expr[valid, , drop = FALSE]
      genes_valid <- gene_symbols[valid]
      expr <- avereps(expr_valid, ID = genes_valid)
    }
  }

  expr_list[[platform]] <- expr
  pdata$platform <- platform
  clin_list[[platform]] <- pdata
}

# 合并同平台（若多个）
# 简单合并：不同平台保持分开保存，并分别输出；并额外输出一个合并后的（以交集基因为主）

# 保存各平台单独结果
for (nm in names(expr_list)) {
  saveRDS(expr_list[[nm]], file = file.path(out_proc_dir, sprintf("GSE14520_expr_%s.rds", nm)))
}
for (nm in names(clin_list)) {
  saveRDS(clin_list[[nm]], file = file.path(out_proc_dir, sprintf("GSE14520_clinical_%s.rds", nm)))
}

# 试图挑选样本数最多的平台作为主平台
platform_sizes <- sapply(expr_list, function(x) ncol(x))
main_platform <- names(which.max(platform_sizes))
expr_main <- expr_list[[main_platform]]
clin_main <- clin_list[[main_platform]]

saveRDS(expr_main, file = file.path(out_proc_dir, "GSE14520_expr.rds"))
saveRDS(clin_main, file = file.path(out_proc_dir, "GSE14520_clinical.rds"))

# log2 转换（若未 log）
q <- quantile(expr_main, probs = 0.99, na.rm = TRUE)
if (q > 50) {
  expr_log2 <- log2(expr_main + 1)
  saveRDS(expr_log2, file = file.path(out_proc_dir, "GSE14520_expr_log2.rds"))
  message("[GSE14520] 已进行 log2 转换并保存为 GSE14520_expr_log2.rds")
} else {
  message("[GSE14520] 看起来已是对数尺度，跳过 log2 转换")
}

# 打印临床字段，帮助人工确认生存与分期信息
cat("\n[临床字段预览 - 前30列] \n")
clin_cols <- colnames(clin_main)
print(head(clin_cols, 30))

# 尝试在常见字段中提取关键信息
survival_keywords <- c("survival", "overall", "os", "dfs", "ttr", "time", "event", "status", "recurrence")
stage_keywords <- c("stage", "grade", "tumor", "tnm")

found_surv <- clin_cols[str_detect(tolower(clin_cols), paste(survival_keywords, collapse = "|"))]
found_stage <- clin_cols[str_detect(tolower(clin_cols), paste(stage_keywords, collapse = "|"))]

cat("\n[可能的生存相关字段]\n")
print(found_surv)
cat("\n[可能的分期/分级相关字段]\n")
print(found_stage)

cat("\n[GSE14520] 处理完成。主平台: ", main_platform, "\n", sep = "")
