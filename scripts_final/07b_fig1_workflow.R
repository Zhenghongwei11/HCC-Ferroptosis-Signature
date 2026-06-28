#!/usr/bin/env Rscript

# 07b_fig1_workflow.R
# Rebuild Figure 1 as a publication-ready study-design workflow.

if (!interactive()) pdf(NULL)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

detect_base_dir <- function() {
  candidates <- c(".", "HCC_Ferroptosis_Project")
  for (d in candidates) {
    if (dir.exists(file.path(d, "results")) && dir.exists(file.path(d, "plots"))) {
      return(d)
    }
  }
  stop("Cannot locate project root.")
}

base_dir <- detect_base_dir()
pub_dir <- file.path(base_dir, "plots", "publication")
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

clr_bg <- "#FBFCFE"
clr_header <- "#1F4E79"
clr_body <- "#FFFFFF"
clr_border <- "#AEBFCB"
clr_text <- "#24323D"
clr_arrow <- "#4A647A"

modules <- data.frame(
  id = c("M1", "M2", "M3", "M4", "M5", "M6"),
  x = c(0.24, 0.76, 0.50, 0.24, 0.76, 0.50),
  y = c(0.84, 0.84, 0.60, 0.35, 0.35, 0.11),
  w = c(0.44, 0.44, 0.94, 0.44, 0.44, 0.94),
  h = c(0.205, 0.205, 0.175, 0.215, 0.215, 0.135),
  header_h = c(0.045, 0.045, 0.042, 0.045, 0.045, 0.038),
  title = c(
    "1. Cohorts and Reference Data",
    "2. Candidate Gene Derivation",
    "3. Prognostic Signature Training",
    "4. External Validation",
    "5. Prediction and Biological Context",
    "6. Interpretation and Next Steps"
  ),
  body = c(
    "Discovery cohort: GSE14520 (n=221; OS)\nExternal cohorts: TCGA-LIHC, ICGC-LIRI-JP,\nGEO (GSE76427/GSE10143/GSE27150)\nFerroptosis reference: FerrDb V2",
    "Differential expression in GSE14520\n(limma Bioconductor package)\nFerrDb genes filtered by DEG significance\nFinal modeling candidate set: 125 genes",
    "Univariate Cox and LASSO-Cox (10-fold CV)\nNine-gene fixed-coefficient risk score\nInternal evaluation: KM survival,\ntime-dependent ROC (1/3/5 years), bootstrap C-index",
    "Testing in TCGA-LIHC, ICGC-LIRI-JP,\nand GEO cohorts\nMetrics: C-index, time-dependent AUC,\n3-year calibration, IPCW Brier\nRobustness: PH checks, RMST,\nrandom-signature benchmark",
    "Prediction-model context: multivariable Cox,\nnomogram, decision-curve analysis\nBiological context: immune signatures,\ncheckpoints, TCGA mutation/CNV/subtypes",
    "Findings support further evaluation\nof HCC risk stratification.\nProspective validation with predefined\nthresholds is required."
  ),
  body_size = c(3.9, 3.8, 3.95, 3.45, 3.55, 3.65),
  stringsAsFactors = FALSE
)

modules$xmin <- modules$x - modules$w / 2
modules$xmax <- modules$x + modules$w / 2
modules$ymin <- modules$y - modules$h / 2
modules$ymax <- modules$y + modules$h / 2
modules$body_ymax <- modules$ymax - modules$header_h
modules$header_y <- (modules$body_ymax + modules$ymax) / 2
modules$body_y <- (modules$ymin + modules$body_ymax) / 2

arrows <- data.frame(
  x = c(0.43, 0.70, 0.40, 0.60, 0.33, 0.67),
  y = c(0.84, 0.73, 0.50, 0.50, 0.25, 0.25),
  xend = c(0.57, 0.55, 0.28, 0.72, 0.44, 0.56),
  yend = c(0.84, 0.67, 0.46, 0.46, 0.17, 0.17)
)

p <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = clr_bg, color = NA) +
  geom_rect(
    data = modules,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = body_ymax),
    fill = clr_body, color = clr_border, linewidth = 0.7
  ) +
  geom_rect(
    data = modules,
    aes(xmin = xmin, xmax = xmax, ymin = body_ymax, ymax = ymax),
    fill = clr_header, color = clr_header, linewidth = 0.7
  ) +
  geom_text(
    data = modules,
    aes(x = x, y = header_y, label = title),
    color = "white", fontface = "bold", size = 4.55, lineheight = 1.0
  ) +
  geom_text(
    data = modules,
    aes(x = x, y = body_y, label = body, size = body_size),
    color = clr_text, lineheight = 1.11, show.legend = FALSE
  ) +
  scale_size_identity() +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = clr_arrow, linewidth = 0.9,
    arrow = arrow(length = unit(0.22, "cm"), type = "closed")
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 16, 8, 16))

out_pdf <- file.path(pub_dir, "Figure1_workflow.pdf")
out_png <- file.path(pub_dir, "Figure1_workflow.png")

ggsave(out_pdf, p, width = 13.5, height = 8.5, device = cairo_pdf, bg = "white")
ggsave(out_png, p, width = 13.5, height = 8.5, dpi = 350, bg = "white")

message("Figure 1 workflow updated: ", out_png)
