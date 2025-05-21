#!/usr/bin/env Rscript
##############################################################################
## Author:  Zhong Guorui
## Created: 2025-05-20
##############################################################################

# Load required libraries
library(optparse)

# Parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Input directory containing SEG files", metavar = "character"
    ),
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output file path for combined data", metavar = "character"
    ),
    make_option(c("-p", "--pattern"),
        type = "character", default = "\\.seg$",
        help = "Pattern to match SEG files (default: \".seg\")", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    stop("Input directory must be specified with -i or --input")
}

# Function to collect and combine SEG files
collect_seg_files <- function(
    input_dir, output_file = NULL, pattern = "\\.seg$") {
    # Find all SEG files
    seg_files <- list.files(
        input_dir,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(seg_files) == 0) {
        stop("No SEG files found in ", input_dir, " matching pattern '", pattern, "'")
    }

    message("Found ", length(seg_files), " SEG files")

    # Default output file name
    if (is.null(output_file)) {
        output_file <- file.path(input_dir, "cohort_combined.seg")
    }

    # Read and combine all files
    combined_df <- data.frame()
    failed_files <- character()

    for (i in seq_along(seg_files)) {
        file <- seg_files[i]
        message("Processing (", i, "/", length(seg_files), "): ", basename(file))

        # Try to read file
        tryCatch(
            {
                df <- read.table(
                    file,
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )

                # Combine with existing data
                if (nrow(combined_df) == 0) {
                    combined_df <- df
                } else {
                    # Determine common columns
                    common_cols <- intersect(
                        colnames(combined_df),
                        colnames(df)
                    )

                    if (length(common_cols) == 0) {
                        warning(
                            "No common columns found with file: ",
                            basename(file)
                        )

                        failed_files <- c(failed_files, basename(file))

                        next
                    }

                    combined_df <- rbind(
                        combined_df[, common_cols],
                        df[, common_cols]
                    )
                }
            },
            error = function(e) {
                warning(
                    "Failed to process file: ", basename(file), " - ",
                    e$message
                )

                failed_files <<- c(failed_files, basename(file))
            }
        )
    }

    if (length(failed_files) > 0) {
        message("Failed to process ", length(failed_files), " files: ", paste(failed_files, collapse = ", "))
    }

    if (nrow(combined_df) == 0) {
        stop("No data could be combined from the input files")
    }

    # Write combined file
    write.table(
        combined_df,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    message(
        "Saved combined file with ", nrow(combined_df),
        " rows and ", ncol(combined_df), " columns: ", output_file
    )

    # # Create filtered version (non-neutral segments)
    # if ("SVTYPE" %in% colnames(combined_df)) {
    #     filtered_output_file <- sub("\\.seg$", "_filtered.seg", output_file)
    #     filtered_df <- combined_df[combined_df$SVTYPE != "NEUTR", ]

    #     write.table(
    #         filtered_df,
    #         file = filtered_output_file,
    #         sep = "\t",
    #         quote = FALSE,
    #         row.names = FALSE
    #     )

    #     message(
    #         "Saved filtered combined file with ", nrow(filtered_df),
    #         " rows and ", ncol(filtered_df), " columns: ", filtered_output_file
    #     )

    #     return(list(
    #         complete = output_file,
    #         filtered = filtered_output_file
    #     ))
    # }

    return(output_file)
}

# If script is being run directly (not sourced)
if (!interactive()) {
    input_dir <- opt$input
    output_file <- opt$output
    pattern <- opt$pattern

    collected_files <- collect_seg_files(input_dir, output_file, pattern)

    if (is.list(collected_files)) {
        message("Collection complete!")
        message("Complete file: ", collected_files$complete)
        message("Filtered file: ", collected_files$filtered)
    } else {
        message("Collection complete! Final file: ", collected_files)
    }
}

# For testing, comment this out when not testing
input_dir <- "/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets"
collect_seg_files(input_dir)
