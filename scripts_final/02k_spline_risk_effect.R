#!/usr/bin/env Rscript

# 02k_spline_risk_effect.R
# Continuous risk-score effect (per SD) + restricted cubic spline visualization.
#
# Outputs:
# - results/risk_score_continuous_effect.csv
# - plots/supplementary/FigureS_risk_score_spline.pdf
# - plots/supplementary/FigureS_risk_score_spline.png
#
# Usage:
#   Rscript scripts_final/02k_spline_risk_effect.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(rms)
  library(ggplot2)
})

res_dir <- "results"
plot_dir <- "plots/supplementary"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

read_risk <- function(path, dataset, time_col = NULL, status_col = "status") {
  if (!file.exists(path)) stop(sprintf("Missing file: %s", path))
  df0 <- read.csv(path, check.names = FALSE)

  if (is.null(time_col)) {
    time_col <- if ("time_months" %in% colnames(df0)) "time_months" else if ("time" %in% colnames(df0)) "time" else NULL
  }
  if (is.null(time_col)) stop(sprintf("No time column found in %s", path))
  if (!all(c("risk_score", time_col, status_col) %in% colnames(df0))) {
    stop(sprintf("Required columns missing in %s", path))
  }

  df <- df0 %>%
    transmute(
      dataset = dataset,
      time_months = as.numeric(.data[[time_col]]),
      status = as.integer(.data[[status_col]]),
      risk_score = as.numeric(risk_score)
    ) %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))

  df
}

extract_per_sd <- function(df) {
  z <- as.numeric(scale(df$risk_score))
  fit <- coxph(Surv(time_months, status) ~ z, data = df)
  s <- summary(fit)
  beta <- as.numeric(coef(fit)[1])
  se <- as.numeric(s$coef[1, "se(coef)"])
  hr <- exp(beta)
  ci_l <- exp(beta - 1.96 * se)
  ci_u <- exp(beta + 1.96 * se)
  p <- as.numeric(s$coef[1, "Pr(>|z|)"])
  tibble(
    Dataset = unique(df$dataset),
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    HR_per_1SD = hr,
    CI95_L = ci_l,
    CI95_U = ci_u,
    P_value = p
  )
}

fit_spline <- function(df) {
  dd <- datadist(df)
  old <- options(datadist = "dd")
  old_dd <- NULL
  if (exists("dd", envir = .GlobalEnv, inherits = FALSE)) {
    old_dd <- get("dd", envir = .GlobalEnv, inherits = FALSE)
  }
  assign("dd", dd, envir = .GlobalEnv)
  on.exit({
    options(old)
    if (is.null(old_dd)) {
      rm("dd", envir = .GlobalEnv)
    } else {
      assign("dd", old_dd, envir = .GlobalEnv)
    }
  }, add = TRUE)

  fit <- cph(Surv(time_months, status) ~ rcs(risk_score, 4), data = df, x = TRUE, y = TRUE, surv = TRUE)
  p <- as.data.frame(Predict(fit, risk_score, ref.zero = TRUE, fun = exp, conf.int = 0.95))
  p <- p %>%
    transmute(
      dataset = unique(df$dataset),
      risk_score = risk_score,
      HR = yhat,
      CI95_L = lower,
      CI95_U = upper
    )
  list(fit = fit, pred = p)
}

train <- read_risk(file.path(res_dir, "risk_score_data.csv"), dataset = "GSE14520 (Discovery)", time_col = "time")
tcga <- read_risk(file.path(res_dir, "TCGA_LIHC_risk_score.csv"), dataset = "TCGA-LIHC", time_col = "time_months")

effect_tab <- bind_rows(
  extract_per_sd(train),
  extract_per_sd(tcga)
)
write.csv(effect_tab, file.path(res_dir, "risk_score_continuous_effect.csv"), row.names = FALSE)
message("✅ Wrote: results/risk_score_continuous_effect.csv")

s1 <- fit_spline(train)$pred
s2 <- fit_spline(tcga)$pred
pred <- bind_rows(s1, s2)

pred$dataset <- factor(pred$dataset, levels = c("GSE14520 (Discovery)", "TCGA-LIHC"))

p <- ggplot(pred, aes(x = risk_score, y = HR)) +
  geom_ribbon(aes(ymin = CI95_L, ymax = CI95_U), fill = "#6FA8B5", alpha = 0.25) +
  geom_line(color = "#2F4F5F", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#9FB1B6", linewidth = 0.5) +
  facet_wrap(~dataset, scales = "free_x") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Risk score",
    y = "Hazard ratio (ref = median)",
    title = "Continuous association between risk score and overall survival"
  )

pdf(file.path(plot_dir, "FigureS_risk_score_spline.pdf"), width = 10, height = 5)
print(p)
dev.off()
message("✅ Wrote: plots/supplementary/FigureS_risk_score_spline.pdf")

ggsave(file.path(plot_dir, "FigureS_risk_score_spline.png"), p, width = 10, height = 5, dpi = 300, bg = "white")
message("✅ Wrote: plots/supplementary/FigureS_risk_score_spline.png")
