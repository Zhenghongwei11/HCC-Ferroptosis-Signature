# 00_setup_env.R
# Optional package bootstrap for local R environments.

options(device = pdf)
graphics.off()
Sys.setenv(DISPLAY = "")

cran_repo <- Sys.getenv("CRAN_REPO", unset = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(repos = c(CRAN = cran_repo))

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

base_pkgs <- c(
  "tidyverse", "ggplot2", "pheatmap", "ggpubr", "RColorBrewer",
  "VennDiagram", "cowplot", "gridExtra", "data.table",
  "GEOquery", "limma", "survival", "survminer", "timeROC",
  "AnnotationDbi", "hgu133a2.db", "glmnet", "pROC", "rms",
  "UCSCXenaTools"
)

for (pkg in base_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    message("Installing: ", pkg)
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }, error = function(e) {
      message("[WARN] Failed to install ", pkg, ": ", e$message)
    })
  }
}

message("Package bootstrap completed.")
