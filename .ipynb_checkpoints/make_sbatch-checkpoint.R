#!/usr/bin/env Rscript

## Minimal dependencies: only yaml
suppressPackageStartupMessages({
  library(yaml)
})

# -----------------------------
# Helper: simple CLI arg parser
# -----------------------------
parse_args <- function(args) {

  # ---- Default values for ALL fields ----
  res <- list(
    config  = "config.yaml",          # default sample list
    submit   = FALSE                    # default: do NOT submit
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    # Flag without value
    if (key == "--submit") {
      res$submit <- TRUE
      i <- i + 1
      next
    }
    # Key-value argument
    if (i == length(args)) {
      stop("Argument ", key, " is missing a value.")
    }
    value <- args[i + 1]
    if (key == "--config") {
      res$config <- value
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 2
  }
  return(res)
}

# --------------------------------
# Helper: template filling function
# --------------------------------
fill_template <- function(template_string, vars) {
  out <- template_string
  for (nm in names(vars)) {
    placeholder <- paste0("{{", nm, "}}")
    out <- gsub(placeholder, vars[[nm]], out, fixed = TRUE)
  }
  return(out)
}

# ---------------
# Main script body
# ---------------
args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(args)

message("Config:  ", opts$config)
message("Submit:         ", opts$submit)

# Create main output folder if not exist
folders = c("./fastp", "./star", "./featurecount", "./deseq2")
for (f in folders){
    if (!dir.exists(f)) {
        dir.create(f)
    }
}

# Read configs
cfg <- read_yaml(opts$config)

if (is.null(cfg$samples) || length(cfg$samples) == 0) {
  stop("No samples found under 'samples' key in ", opts$samples)
}

# Read star sbatch template as a single string 
template_text <- readLines("star_template.sbatch")
template_text <- paste(template_text, collapse = "\n")
# Build a common system-level variable list
# (Guard with `[[` and %||% if you want more robustness)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Set index
if (cfg$genome == "hg38") {
    star_index = cfg$hg38_index
    gtf = cfg$hg38_gtf
} else if (cfg$genome == "mm10") {
    star_index = cfg$mm10_index
    gtf = cfg$mm10_gtf
} else {
    stop("Choose genome from hg38/mm10. Unsupported genome: ", cfg$genome)
}

# Check if sample ids are unique
if (length(unique(sapply(cfg$samples, function(x) x$sample_id))) != length(cfg$samples)) {
  stop("Sample IDs must be unique.")
}

# Loop over samples
samples_list <- cfg$samples
generated_files <- character(0)
job_ids <- character(0)     # <--- store sbatch job IDs
bam_files <- character(0)   # <--- store BAM paths for featureCounts
                         
for (s in samples_list) {
  sample_id <- s$sample_id
  if (is.null(sample_id)) next
  # Combine system-level and sample-level variables
  vars <- list(
      SAMPLE_ID = sample_id,
      FASTQ1    = s$fastq1 %||% "",
      FASTQ2    = s$fastq2 %||% "",
      STAR_PATH = cfg$star_path,
      STAR_INDEX = star_index
    )
  rendered <- fill_template(template_text, vars)
  sbatch_file <- file.path(paste0("STEP01_", sample_id, ".sh"))
  writeLines(rendered, sbatch_file)
  message("Generated: ", sbatch_file)
  generated_files <- c(generated_files, sbatch_file)
  # Construct BAM path used later for featureCounts
  bam_path <- file.path("./star/", paste0(sample_id, "_Aligned.sortedByCoord.out.bam"))
  bam_files <- c(bam_files, bam_path)

  if (opts$submit) {
    cmd <- paste("sbatch", sbatch_file)
    message("Submitting: ", cmd)
    sbatch_out <- system(cmd, intern = TRUE)
    message("  -> ", sbatch_out)
      # Parse job ID from "Submitted batch job 123456"
    jid <- sub("Submitted batch job[[:space:]]+([0-9]+).*", "\\1", sbatch_out)
    job_ids <- c(job_ids, jid)
  } else {
    message("Not submitting (use --submit to enable).")
  }
}

message("Generated ", length(generated_files), " sbatch file(s).")

# -----------------------------
# STEP02: featureCounts sbatch
# -----------------------------
# Only generate if we actually have BAMs
if (length(bam_files) > 0) {
  template_text <- readLines("featurecounts_template.sbatch")
  template_text <- paste(template_text, collapse = "\n")
  vars <- list(
      FEATURECOUNTS_PATH = cfg$featurecounts_path,
      GTF = gtf,
      BAMS = paste(bam_files, collapse = " ")
  )
  rendered <- fill_template(template_text, vars)
  fc_file <- file.path(paste0("STEP02_featurecounts.sh"))
  writeLines(rendered, fc_file)
  message("Generated featureCounts sbatch: ", sbatch_file)
  if (opts$submit) {
    dep <- paste(job_ids, collapse = ":")
    dep_flag <- paste0("--dependency=afterok:", dep)
    fc_cmd <- paste("sbatch", dep_flag, fc_file)
    message("Submitting featureCounts job: ", fc_cmd)
    system(fc_cmd)
  }
}

message("Done. Generated featureCounts files!")
