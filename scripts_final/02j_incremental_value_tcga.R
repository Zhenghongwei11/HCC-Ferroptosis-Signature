#!/usr/bin/env Rscript

# 02j_incremental_value_tcga.R
# Incremental value of the gene-based risk score beyond TCGA clinical covariates.
#
# Outputs:
# - results/tcga_incremental_value.csv
#
# Usage:
#   Rscript scripts_final/02j_incremental_value_tcga.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

tcga_path <- file.path(res_dir, "TCGA_LIHC_risk_score.csv")
if (!file.exists(tcga_path)) stop("Missing TCGA risk-score file: results/TCGA_LIHC_risk_score.csv")

df0 <- read.csv(tcga_path, check.names = FALSE)
needed <- c("risk_score", "time_months", "status", "age_at_diagnosis", "gender", "ajcc_pathologic_stage", "tumor_grade")
missing <- setdiff(needed, colnames(df0))
if (length(missing) > 0) {
  stop(paste0("Missing columns in TCGA file: ", paste(missing, collapse = ", ")))
}

parse_stage <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- ifelse(is.na(x) | x == "", NA_character_, x)
  x <- gsub("^Stage\\s+", "", x, ignore.case = TRUE)
  x <- gsub("[A-Za-z].*$", "", x) # IIIA -> III
  x <- trimws(x)
  x
}

df <- df0 %>%
  transmute(
    time_months = as.numeric(time_months),
    status = as.integer(status),
    risk_score = as.numeric(risk_score),
    age_years = as.numeric(age_at_diagnosis) / 365.25,
    gender = as.factor(as.character(gender)),
    stage = parse_stage(ajcc_pathologic_stage),
    grade = as.character(tumor_grade)
  ) %>%
  filter(!is.na(time_months), !is.na(status), !is.na(risk_score)) %>%
  mutate(
    gender = forcats::fct_na_value_to_level(gender, level = "Unknown"),
    stage = factor(stage, levels = c("I", "II", "III", "IV")),
    stage = forcats::fct_na_value_to_level(stage, level = "Unknown"),
    grade = factor(grade, levels = c("G1", "G2", "G3", "G4")),
    grade = forcats::fct_na_value_to_level(grade, level = "Unknown")
  )

df <- df %>% filter(!is.na(age_years) & is.finite(age_years))

surv_obj <- Surv(df$time_months, df$status)

cox_gene <- coxph(surv_obj ~ risk_score, data = df)
cox_clin <- coxph(surv_obj ~ age_years + gender + stage + grade, data = df)
cox_comb <- coxph(surv_obj ~ risk_score + age_years + gender + stage + grade, data = df)

extract_metrics <- function(fit) {
  s <- summary(fit)
  cidx <- as.numeric(s$concordance[1])
  cse <- as.numeric(s$concordance[2])
  tibble(
    C_index = cidx,
    C_index_SE = cse,
    LogLik = as.numeric(logLik(fit)[1]),
    AIC = as.numeric(AIC(fit))
  )
}

tab <- bind_rows(
  tibble(Model = "Gene-only") %>% bind_cols(extract_metrics(cox_gene)),
  tibble(Model = "Clinical-only") %>% bind_cols(extract_metrics(cox_clin)),
  tibble(Model = "Combined") %>% bind_cols(extract_metrics(cox_comb))
) %>%
  mutate(
    Dataset = "TCGA-LIHC",
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE)
  ) %>%
  select(Dataset, Sample_N, Events, Model, C_index, C_index_SE, LogLik, AIC)

lrt_gene_vs_comb <- anova(cox_gene, cox_comb, test = "LRT")
lrt_clin_vs_comb <- anova(cox_clin, cox_comb, test = "LRT")

lrt_row <- tibble(
  Dataset = "TCGA-LIHC",
  Contrast = c("Gene-only vs Combined", "Clinical-only vs Combined"),
  LRT_df = c(lrt_gene_vs_comb$Df[2], lrt_clin_vs_comb$Df[2]),
  LRT_chisq = c(lrt_gene_vs_comb$Chisq[2], lrt_clin_vs_comb$Chisq[2])
) %>%
  mutate(
    LRT_p = stats::pchisq(LRT_chisq, df = LRT_df, lower.tail = FALSE)
)

write.csv(tab, file.path(res_dir, "tcga_incremental_value.csv"), row.names = FALSE)
write.csv(lrt_row, file.path(res_dir, "tcga_incremental_value_lrt.csv"), row.names = FALSE)

extract_hr <- function(fit, term) {
  s <- summary(fit)
  if (!(term %in% rownames(s$coef))) stop(sprintf("Term not found in model: %s", term))
  beta <- as.numeric(s$coef[term, "coef"])
  se <- as.numeric(s$coef[term, "se(coef)"])
  hr <- exp(beta)
  ci_l <- exp(beta - 1.96 * se)
  ci_u <- exp(beta + 1.96 * se)
  p <- as.numeric(s$coef[term, "Pr(>|z|)"])
  tibble(
    Term = term,
    HR = hr,
    CI95_L = ci_l,
    CI95_U = ci_u,
    P_value = p
  )
}

sd_risk <- sd(df$risk_score, na.rm = TRUE)
adj_risk <- extract_hr(cox_comb, "risk_score") %>%
  mutate(
    Dataset = "TCGA-LIHC",
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    RiskScore_SD = sd_risk,
    HR_per_1SD = exp(log(HR) * sd_risk),
    CI95_L_per_1SD = exp(log(CI95_L) * sd_risk),
    CI95_U_per_1SD = exp(log(CI95_U) * sd_risk)
  ) %>%
  select(
    Dataset, Sample_N, Events,
    RiskScore_SD, HR, CI95_L, CI95_U, P_value,
    HR_per_1SD, CI95_L_per_1SD, CI95_U_per_1SD
  )

write.csv(adj_risk, file.path(res_dir, "tcga_combined_risk_score_adjusted_hr.csv"), row.names = FALSE)
message("✅ Wrote: results/tcga_combined_risk_score_adjusted_hr.csv")

message("✅ Wrote: results/tcga_incremental_value.csv")
message("✅ Wrote: results/tcga_incremental_value_lrt.csv")
