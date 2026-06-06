#!/usr/bin/env Rscript

# run_complete_pipeline.R
# HCC ferroptosis signature main pipeline for the current submission.
# Usage: Rscript scripts_final/run_complete_pipeline.R

# ============================================================================
# 禁用交互式图形设备
# ============================================================================
options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

message("================================================================================")
message("HCC ferroptosis signature analysis - current final pipeline")
message("================================================================================")
message("开始时间: ", Sys.time())

# 检查工作目录
if (!dir.exists("data/processed")) {
  stop("错误: 请在项目根目录运行此脚本")
}

# ============================================================================
# Stage 1: 多队列整合与验证
# ============================================================================
message("\n[Stage 1] 数据准备与基础分析...")
source("scripts_final/00_setup_env.R")
# 准备 HCC 特异 ferroptosis 基因集
source("scripts_final/00c_prepare_hcc_ferro_genes.R")
source("scripts_final/01b_download_GSE14520.R")
source("scripts_final/02b_multi_cohort_DEG.R")

# ============================================================================
# Stage 2: 预后模型与临床效用
# ============================================================================
message("\n[Stage 2] 预后模型构建与验证...")
source("scripts_final/02c_prognostic_model.R")
source("scripts_final/02d_nomogram_calibration_dca.R")
source("scripts_final/02e_external_validation.R")

# Diagnostics/robustness for final manuscript
source("scripts_final/02h_model_diagnostics.R")

# Calibration, prediction error, incremental value, and spline visualization
source("scripts_final/02i_external_calibration_brier.R")
source("scripts_final/02j_incremental_value_tcga.R")
source("scripts_final/02k_spline_risk_effect.R")
source("scripts_final/02m_tcga_multiomics_characterization.R")

# ============================================================================
# Stage 3: 免疫微环境分析
# ============================================================================
message("\n[Stage 3] 免疫微环境分析...")
source("scripts_final/03a_immune_infiltration.R")
source("scripts_final/03b_immune_checkpoint.R")

# ============================================================================
# Stage 4: 最终出版级图表生成
# ============================================================================
message("\n[Stage 4] 生成最终出版级图表 (Figures 2-5)...")
source("scripts_final/07_final_publication_figures.R")
source("scripts_final/07b_fig1_workflow.R")

message("\n================================================================================")
message("✅ 分析管道全部完成！")
message("请查看 'plots/publication/' 目录获取最终图表。")
message("================================================================================")
