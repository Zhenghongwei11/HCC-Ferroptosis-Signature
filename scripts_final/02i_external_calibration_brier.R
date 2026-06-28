#!/usr/bin/env Rscript

# 02i_external_calibration_brier.R
# External calibration + Brier score evaluation
#
# - Uses the discovery cohort to estimate baseline cumulative hazard for the fixed linear predictor
#   (risk_score as Cox linear predictor; coefficient fixed to 1 via offset).
# - Applies the discovery baseline hazard to all cohorts to generate predicted risk at 1/3/5 years
#   and evaluates Brier score with IPCW (pec::sbrier).
# - Produces 3-year calibration plots (pec::calPlot) per cohort.
#
# Usage: Rscript scripts_final/02i_external_calibration_brier.R

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(pec)
})

res_dir <- "results"
plot_dir <- "plots/supplementary"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

read_risk_csv <- function(path, dataset) {
  if (!file.exists(path)) stop(sprintf("Missing risk score file for %s: %s", dataset, path))
  df <- read.csv(path, check.names = FALSE)

  if ("time" %in% colnames(df) && !"time_months" %in% colnames(df)) {
    df$time_months <- df$time
  }
  if (!all(c("risk_score", "time_months", "status") %in% colnames(df))) {
    stop(sprintf("Required columns not found in %s: need risk_score,time_months,status", path))
  }

  df <- df %>%
    transmute(
      dataset = dataset,
      time_months = as.numeric(time_months),
      status = as.integer(status),
      risk_score = as.numeric(risk_score)
    ) %>%
    filter(!is.na(time_months), !is.na(status), !is.na(risk_score))

  df
}

# 1) Fit offset Cox in discovery to estimate baseline cumulative hazard for the fixed linear predictor.
train_path <- file.path(res_dir, "risk_score_data.csv")
train_df <- read_risk_csv(train_path, dataset = "GSE14520 (Discovery)")
train_surv <- Surv(train_df$time_months, train_df$status)

cox_offset <- coxph(train_surv ~ offset(risk_score), data = train_df, x = TRUE, y = TRUE)
bh <- basehaz(cox_offset, centered = FALSE) %>% as_tibble()
if (!all(c("time", "hazard") %in% colnames(bh))) stop("Unexpected basehaz() output.")

H0_step <- stats::approxfun(
  x = bh$time,
  y = bh$hazard,
  method = "constant",
  f = 0,
  rule = 2
)

predict_surv_prob <- function(lp, t_months) {
  H0t <- as.numeric(H0_step(t_months))
  exp(-H0t * exp(lp))
}

ipcw_brier <- function(time_months, status, pred_risk, t_months) {
  # IPCW Brier score for right-censored survival data at a fixed horizon t.
  # status: 1=event, 0=censored
  # pred_risk: predicted event probability by time t (0..1)
  if (length(time_months) != length(status) || length(status) != length(pred_risk)) {
    stop("ipcw_brier: input length mismatch")
  }
  df <- tibble(
    time_months = as.numeric(time_months),
    status = as.integer(status),
    pred_risk = as.numeric(pred_risk)
  ) %>%
    filter(!is.na(time_months), !is.na(status), !is.na(pred_risk))

  # Censoring distribution G(t) = P(C >= t)
  # Treat censoring as event: 1-status
  gfit <- survfit(Surv(df$time_months, 1 - df$status) ~ 1)
  gstep <- stats::approxfun(
    x = gfit$time,
    y = gfit$surv,
    method = "constant",
    f = 0,
    rule = 2
  )

  tt <- t_months
  y <- as.integer(df$time_months <= tt & df$status == 1)
  # weights:
  # - event before tt: 1/G(time)
  # - event-free at tt: 1/G(tt)
  # - censored before tt: 0
  w <- rep(0, nrow(df))
  is_event_before <- df$time_months <= tt & df$status == 1
  is_at_risk_tt <- df$time_months > tt

  g_event <- pmax(as.numeric(gstep(df$time_months[is_event_before])), 1e-12)
  g_tt <- pmax(as.numeric(gstep(tt)), 1e-12)

  w[is_event_before] <- 1 / g_event
  w[is_at_risk_tt] <- 1 / g_tt

  mean(w * (y - df$pred_risk) ^ 2, na.rm = TRUE)
}

cohorts <- list(
  "TCGA-LIHC" = file.path(res_dir, "TCGA_LIHC_risk_score.csv"),
  "GSE76427" = file.path(res_dir, "GSE76427_risk_score.csv"),
  "GSE10143-HCC" = file.path(res_dir, "GSE10143_HCC_risk_score.csv"),
  "GSE27150" = file.path(res_dir, "GSE27150_risk_score.csv"),
  "ICGC-LIRI-JP (HCCDB18)" = file.path(res_dir, "HCCDB18_LIRIJP_risk_score.csv")
)

times <- c(12, 36, 60)

# 2) Brier scores (IPCW) at 1/3/5 years per cohort.
brier_rows <- list()
for (nm in names(cohorts)) {
  df <- read_risk_csv(cohorts[[nm]], dataset = nm)
  surv_obj <- Surv(df$time_months, df$status)

  max_t <- max(df$time_months, na.rm = TRUE)
  eval_times <- times[times < max_t]

  out <- tibble(
    Dataset = nm,
    Sample_N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    Brier_1y = NA_real_,
    Brier_3y = NA_real_,
    Brier_5y = NA_real_
  )

  for (t in eval_times) {
    pred_surv <- predict_surv_prob(df$risk_score, t)
    pred_risk <- 1 - pred_surv
    b <- ipcw_brier(df$time_months, df$status, pred_risk, t)
    if (t == 12) out$Brier_1y <- b
    if (t == 36) out$Brier_3y <- b
    if (t == 60) out$Brier_5y <- b
  }

  brier_rows[[nm]] <- out
}

brier_df <- bind_rows(brier_rows)
write.csv(brier_df, file.path(res_dir, "external_validation_brier.csv"), row.names = FALSE)
message("Wrote: results/external_validation_brier.csv")

# 3) 3-year calibration plots (risk at 36 months) per cohort.
cal_time <- 36
draw_calibration_panels <- function() {
  par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

  for (nm in names(cohorts)) {
    df <- read_risk_csv(cohorts[[nm]], dataset = nm)

    max_t <- max(df$time_months, na.rm = TRUE)
    if (cal_time >= max_t) {
      plot.new()
      title(main = sprintf("%s: 3-year calibration\nInsufficient follow-up", nm), cex.main = 1)
      next
    }

    pred_surv <- predict_surv_prob(df$risk_score, cal_time)
    pred_risk <- 1 - pred_surv

    obj <- list("Predicted risk" = matrix(pred_risk, ncol = 1))
    calPlot(
      object = obj,
      time = cal_time,
      formula = Surv(time_months, status) ~ 1,
      data = df,
      type = "risk",
      pseudo = FALSE,
      q = 10,
      xlab = "Predicted event probability (%)",
      ylab = "Observed event probability (%)",
      legend = FALSE,
      round = TRUE,
      bars = TRUE,
      hanging = TRUE,
      col = "#1F77B4",
      lwd = 2
    )
    title(
      main = sprintf(
        "%s (N=%d, Events=%d)",
        nm, nrow(df), sum(df$status == 1, na.rm = TRUE)
      ),
      cex.main = 1
    )
  }
}

pdf(file.path(plot_dir, "FigureS_external_calibration_3y.pdf"), width = 12, height = 12)
draw_calibration_panels()
dev.off()
message("Wrote: plots/supplementary/FigureS_external_calibration_3y.pdf")

png(file.path(plot_dir, "FigureS_external_calibration_3y.png"), width = 3600, height = 3600, res = 300)
draw_calibration_panels()
dev.off()
message("Wrote: plots/supplementary/FigureS_external_calibration_3y.png")
