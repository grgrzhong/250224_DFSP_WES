#!/usr/bin/env Rscript

# Parse command line arguments
suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
    library(here)
})

# Define command line options
option_list <- list(
    make_option(c("--facet_dir"),
        type = "character", default = NULL,
        help = "Directory containing FACET VCF files", metavar = "PATH"
    ),
    make_option(c("--output"),
        type = "character", default = "facet_cnv_summary.csv",
        help = "Output file path [default: %default]", metavar = "FILE"
    )
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$facet_dir)) {
    stop("--facet_dir argument is required")
}

# Source study library if available (make this optional)
tryCatch(
    {
        source(here::here("lib/R/study_lib.R"))
    },
    error = function(e) {
        message("Note: study_lib.R not found or couldn't be loaded.")
    }
)

# List all VCF files
vcf_files <- list.files(
    path = here("data/cnv/facets"),
    pattern = "\\.vcf\\.gz$",
    recursive = TRUE,
    full.names = TRUE
)

# Function to extract values from VCF header
extract_facet_info <- function(vcf_file) {
    header <- system(paste("zcat", shQuote(vcf_file), "| grep '^##'"), intern = TRUE)
    get_val <- function(key) {
        val <- sub(
            paste0("##", key, "="), "", header[grep(paste0("^##", key, "="), header)]
        )
        if (length(val) == 0) {
            return(NA)
        }
        val
    }
    sample <- basename(dirname(vcf_file))
    tibble(
        sample = sample,
        purity = as.numeric(get_val("purity")),
        ploidy = as.numeric(get_val("ploidy")),
        dipLogR = as.numeric(get_val("dipLogR")),
        est_insert_size = as.numeric(get_val("est_insert_size"))
    )
}

# Extract info for all files
facet_info <- list_rbind(lapply(vcf_files, extract_facet_info))

facet_info <- facet_info |> 
    mutate(
        sample = str_extract(sample, "^[^_vs]+")
    ) |> 
    arrange(sample)

# Save to CSV
write_excel_csv(
    facet_info, 
    file = "facet_cnv_summary.csv"
)

# Print summary
message(paste("Processed", nrow(facet_info), "samples"))
message(paste("Output saved to", opt$output))