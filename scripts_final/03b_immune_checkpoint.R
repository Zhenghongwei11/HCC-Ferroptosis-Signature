#!/usr/bin/env Rscript

# 03b_immune_checkpoint.R
# 免疫检查点分析 + dysfunction/exclusion 相关探索性替代评分
# 使用方法: Rscript scripts_final/03b_immune_checkpoint.R
# 
# 输入: 表达矩阵, Risk Score数据
# 输出: 免疫检查点分析结果, Figure 9

# 禁用交互式图形设备 (防止XQuartz弹出)
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(corrplot)
})

# 设置工作目录
if (!dir.exists("data/processed")) {
  stop("请在项目根目录运行此脚本")
}

proc_dir <- "data/processed"
res_dir  <- "results"
plot_dir <- "plots"

set.seed(123)

message("[免疫检查点] 开始免疫检查点分析...")

# ============================================
# 1. 加载数据
# ============================================
message("[免疫检查点] 加载数据...")

# 使用基因符号版本的表达矩阵
expr_symbol_file <- file.path(proc_dir, "GSE14520_expr_symbol.rds")
if (file.exists(expr_symbol_file)) {
  expr_14 <- readRDS(expr_symbol_file)
  message("[免疫检查点] 使用基因符号版本表达矩阵")
} else {
  stop("请先运行03a_immune_infiltration.R生成基因符号版本表达矩阵")
}

# Risk score data for the final nine-gene model
risk_file <- file.path(res_dir, "risk_score_data.csv")
if (!file.exists(risk_file)) {
  stop(
    "[免疫检查点] 未找到 results/risk_score_data.csv\n",
    "请先运行: Rscript scripts_final/02c_prognostic_model.R\n",
    "注意：论文主线以最终九基因模型为准。"
  )
}
risk_data <- read.csv(risk_file)

message("[免疫检查点] 样本数: ", ncol(expr_14))

# ============================================
# 2. 定义免疫检查点基因
# ============================================
message("[免疫检查点] 定义免疫检查点基因...")

# 免疫检查点基因列表
checkpoint_genes <- list(
  # 抑制性检查点
  "Inhibitory" = c(
    "PDCD1",      # PD-1
    "CD274",      # PD-L1
    "PDCD1LG2",   # PD-L2
    "CTLA4",      # CTLA-4
    "LAG3",       # LAG-3
    "HAVCR2",     # TIM-3
    "TIGIT",      # TIGIT
    "BTLA",       # BTLA
    "VSIR",       # VISTA
    "IDO1",       # IDO1
    "SIGLEC15"    # Siglec-15
  ),
  # 刺激性检查点
  "Stimulatory" = c(
    "CD28",       # CD28
    "ICOS",       # ICOS
    "TNFRSF4",    # OX40
    "TNFRSF9",    # 4-1BB
    "CD40",       # CD40
    "CD40LG",     # CD40L
    "TNFRSF18",   # GITR
    "CD27",       # CD27
    "CD70"        # CD70
  )
)

all_checkpoints <- unlist(checkpoint_genes)

# ============================================
# 3. 提取免疫检查点表达
# ============================================
message("[免疫检查点] 提取免疫检查点表达...")

available_checkpoints <- intersect(all_checkpoints, rownames(expr_14))
message("[免疫检查点] 可用检查点基因: ", length(available_checkpoints), "/", length(all_checkpoints))

checkpoint_expr <- expr_14[available_checkpoints, , drop = FALSE]

# ============================================
# 4. 高低风险组检查点表达差异
# ============================================
message("[免疫检查点] 分析高低风险组差异...")

common_samples <- intersect(risk_data$sample, colnames(checkpoint_expr))
checkpoint_matched <- checkpoint_expr[, common_samples]
risk_matched <- risk_data[match(common_samples, risk_data$sample), ]

diff_results <- data.frame(
  Gene = rownames(checkpoint_matched),
  High_Mean = numeric(nrow(checkpoint_matched)),
  Low_Mean = numeric(nrow(checkpoint_matched)),
  Log2FC = numeric(nrow(checkpoint_matched)),
  P_value = numeric(nrow(checkpoint_matched)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(checkpoint_matched)) {
  high_vals <- as.numeric(checkpoint_matched[i, risk_matched$risk_group == "High"])
  low_vals <- as.numeric(checkpoint_matched[i, risk_matched$risk_group == "Low"])
  
  diff_results$High_Mean[i] <- mean(high_vals, na.rm = TRUE)
  diff_results$Low_Mean[i] <- mean(low_vals, na.rm = TRUE)
  diff_results$Log2FC[i] <- log2((diff_results$High_Mean[i] + 0.01) / (diff_results$Low_Mean[i] + 0.01))
  
  if (length(high_vals) > 3 && length(low_vals) > 3) {
    diff_results$P_value[i] <- wilcox.test(high_vals, low_vals)$p.value
  }
}

diff_results <- diff_results %>% arrange(P_value)
write.csv(diff_results, file.path(res_dir, "checkpoint_risk_diff.csv"), row.names = FALSE)

message("[免疫检查点] 显著差异的检查点: ", sum(diff_results$P_value < 0.05))

# ============================================
# 5. 检查点与Risk Score相关性
# ============================================
message("[免疫检查点] 计算与Risk Score相关性...")

cor_results <- data.frame(
  Gene = rownames(checkpoint_matched),
  Correlation = numeric(nrow(checkpoint_matched)),
  P_value = numeric(nrow(checkpoint_matched))
)

for (i in 1:nrow(checkpoint_matched)) {
  cor_test <- cor.test(as.numeric(checkpoint_matched[i, ]), risk_matched$risk_score,
                       method = "spearman")
  cor_results$Correlation[i] <- cor_test$estimate
  cor_results$P_value[i] <- cor_test$p.value
}

cor_results <- cor_results %>% arrange(P_value)
write.csv(cor_results, file.path(res_dir, "checkpoint_risk_correlation.csv"), row.names = FALSE)

# ============================================
# 6. dysfunction/exclusion 替代评分计算
# ============================================
message("[免疫检查点] 计算dysfunction/exclusion替代评分...")

# 使用功能障碍和排斥相关标志基因的平均表达构建探索性替代指标；
# 并非官方免疫治疗响应预测算法输出或免疫治疗响应验证。

# 计算免疫排斥评分 (基于T细胞功能障碍相关基因)
dysfunction_genes <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
available_dys <- intersect(dysfunction_genes, rownames(expr_14))

if (length(available_dys) > 0) {
  dysfunction_score <- colMeans(expr_14[available_dys, , drop = FALSE], na.rm = TRUE)
} else {
  dysfunction_score <- rep(NA, ncol(expr_14))
}

# 计算免疫排斥评分 (基于CAF和MDSC相关基因)
exclusion_genes <- c("FAP", "ACTA2", "COL1A1", "COL1A2", "TGFB1")
available_exc <- intersect(exclusion_genes, rownames(expr_14))

if (length(available_exc) > 0) {
  exclusion_score <- colMeans(expr_14[available_exc, , drop = FALSE], na.rm = TRUE)
} else {
  exclusion_score <- rep(NA, ncol(expr_14))
}

# dysfunction/exclusion 替代评分 = 标准化功能障碍 + 标准化排斥
surrogate_score <- scale(dysfunction_score) + scale(exclusion_score)
names(surrogate_score) <- colnames(expr_14)

# 保存替代评分
surrogate_data <- data.frame(
  sample = names(surrogate_score),
  dysfunction_score = dysfunction_score,
  exclusion_score = exclusion_score,
  dysfunction_exclusion_score = as.numeric(surrogate_score)
)

# 匹配Risk数据
surrogate_data <- surrogate_data %>%
  left_join(risk_data[, c("sample", "risk_score", "risk_group")], by = "sample")

write.csv(surrogate_data, file.path(res_dir, "dysfunction_exclusion_surrogate.csv"), row.names = FALSE)

# ============================================
# 7. 生成Figure 9: 免疫检查点多面板图
# ============================================
message("[免疫检查点] 生成Figure 9...")

pdf(file.path(plot_dir, "Figure9_immune_checkpoint.pdf"), width = 14, height = 12)

par(mfrow = c(2, 2))

# A: 免疫检查点表达热图
checkpoint_scaled <- t(scale(t(as.matrix(checkpoint_matched))))
checkpoint_scaled[checkpoint_scaled > 2] <- 2
checkpoint_scaled[checkpoint_scaled < -2] <- -2

# 按Risk分组排序
order_idx <- order(risk_matched$risk_group, risk_matched$risk_score)
checkpoint_ordered <- checkpoint_scaled[, order_idx]

# 颜色条
group_colors <- ifelse(risk_matched$risk_group[order_idx] == "High", "red", "blue")

heatmap(checkpoint_ordered,
        main = "A. Immune Checkpoint Expression",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        scale = "none", Colv = NA, cexRow = 0.7,
        ColSideColors = group_colors)

# B: 高低风险组检查点差异
sig_checkpoints <- diff_results %>% filter(P_value < 0.1) %>% head(10)

if (nrow(sig_checkpoints) > 0) {
  diff_mat <- as.matrix(sig_checkpoints[, c("High_Mean", "Low_Mean")])
  rownames(diff_mat) <- sig_checkpoints$Gene
  
  barplot(t(diff_mat), beside = TRUE,
          col = c("red", "blue"), las = 2,
          main = "B. Checkpoint Expression: High vs Low Risk",
          ylab = "Expression Level")
  legend("topright", legend = c("High Risk", "Low Risk"),
         fill = c("red", "blue"))
}

# C: 检查点与Risk Score相关性
top_cor <- head(cor_results, 10)
barplot(top_cor$Correlation, names.arg = top_cor$Gene,
        col = ifelse(top_cor$Correlation > 0, "red", "blue"),
        las = 2, main = "C. Checkpoint - Risk Score Correlation",
        ylab = "Spearman Correlation")
abline(h = 0, lty = 2)

# D: dysfunction/exclusion替代评分与Risk分组
if (!all(is.na(surrogate_data$dysfunction_exclusion_score))) {
  surrogate_complete <- surrogate_data[complete.cases(surrogate_data[, c("dysfunction_exclusion_score", "risk_group")]), ]
  
  boxplot(dysfunction_exclusion_score ~ risk_group, data = surrogate_complete,
          col = c("blue", "red"),
          main = "D. Dysfunction/Exclusion Score by Risk Group",
          xlab = "Risk Group", ylab = "Dysfunction/exclusion score")
  
  # 添加p值
  if (nrow(surrogate_complete) > 10) {
    p_val <- wilcox.test(dysfunction_exclusion_score ~ risk_group, data = surrogate_complete)$p.value
    mtext(paste("p =", format(p_val, digits = 3)), side = 3, line = -2)
  }
}

dev.off()

message("[免疫检查点] Figure 9已保存")

# ============================================
# 8. 基于替代评分的中位数分组摘要
# ============================================
message("[免疫检查点] 汇总dysfunction/exclusion替代评分分组...")

# 基于替代评分中位数构建探索性分组，不代表真实免疫治疗响应。
surrogate_complete <- surrogate_data[complete.cases(surrogate_data$dysfunction_exclusion_score), ]
surrogate_complete$surrogate_group <- ifelse(
  surrogate_complete$dysfunction_exclusion_score < median(surrogate_complete$dysfunction_exclusion_score, na.rm = TRUE),
  "Lower dysfunction/exclusion", "Higher dysfunction/exclusion"
)

# 统计
surrogate_by_risk <- table(surrogate_complete$risk_group, surrogate_complete$surrogate_group)
message("[免疫检查点] 替代评分分组:")
print(surrogate_by_risk)

write.csv(surrogate_by_risk, file.path(res_dir, "dysfunction_exclusion_group_summary.csv"))

# ============================================
# 完成
# ============================================
message("[免疫检查点] 免疫检查点分析完成！")
message("[免疫检查点] 输出文件:")
message("  ", file.path(res_dir, "checkpoint_risk_diff.csv"))
message("  ", file.path(res_dir, "checkpoint_risk_correlation.csv"))
message("  ", file.path(res_dir, "dysfunction_exclusion_surrogate.csv"))
message("  ", file.path(res_dir, "dysfunction_exclusion_group_summary.csv"))
message("  ", file.path(plot_dir, "Figure9_immune_checkpoint.pdf"))
