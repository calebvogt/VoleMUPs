# Render one VoleMUPs R Markdown report to HTML in the project's Synology reports folder.
# usage (from the repo root): Rscript R/render_report.R 01_mup_search
# Modelled on 2027_GenomicMonogamy/R/render_report.R. Needs pandoc: uses RStudio's copy if RSTUDIO_PANDOC is not set.

suppressPackageStartupMessages({library(rmarkdown); library(yaml); library(here)})
args <- commandArgs(trailingOnly = TRUE)
name <- if (length(args)) args[1] else "01_mup_search"
proj <- here()
cfg <- read_yaml(file.path(proj, "config.yml"))
if (file.exists(file.path(proj, "config.local.yml"))) cfg <- modifyList(cfg, read_yaml(file.path(proj, "config.local.yml")))
if (Sys.getenv("RSTUDIO_PANDOC") == "") {
  for (p in c(file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "RStudio", "resources", "app", "bin", "quarto", "bin", "tools"),
              "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools")) {
    if (dir.exists(p)) { Sys.setenv(RSTUDIO_PANDOC = p); break }
  }
}
out_dir <- file.path(dirname(cfg$data_root), "reports")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
render(file.path(proj, "reports", paste0(name, ".Rmd")), output_dir = out_dir, quiet = TRUE,
       intermediates_dir = tempdir(), knit_root_dir = proj)
cat("wrote", file.path(out_dir, paste0(name, ".html")), "\n")
