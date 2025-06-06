## Load required libraries and functions
source(here::here("bin/R/lib/study_lib.R"))

# Define paths
work_dir <- "/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395"
# work_dir <- "/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/NA12878"

mutect2_dirs <- dir_ls(work_dir, pattern = "mutect2_", type = "directory")
mutect2_dirs <- mutect2_dirs[str_detect(basename(mutect2_dirs), "^mutect2_")]

summary_data <- tibble()

for (mutect2_dir in mutect2_dirs) {
    # Extract depth from directory name
    depth <- str_replace(basename(mutect2_dir), "mutect2_", "") |> as.numeric()

    # Define isec directory
    isec_dir <- file.path(mutect2_dir, "isec")

    # Check if isec directory exists
    if (!dir.exists(isec_dir)) {
        warning(paste("Directory does not exist:", isec_dir))
        next
    }

    # Read the summary CSV file
    csv_file <- file.path(isec_dir, "isec_summary.csv")

    if (file.exists(csv_file)) {
        isec_data <- read_csv(csv_file, show_col_types = FALSE)

        # Print the summary for this depth
        cat("=== Summary for Depth", depth, "===\n")
        print(isec_data)
        cat("\n")
    } else {
        warning(paste("CSV file not found for depth", depth, ":", csv_file))
    }

    # Append the summary data to isec_data
    summary_data <- bind_rows(summary_data, isec_data)
}

## Calculate the performance metrics
performance_metrics <- summary_data |>
    mutate(
        # Basic counts
        TP = intersection, # True Positives
        FP = mutect2_unique, # False Positives (called but not in truth)
        FN = truth_unique, # False Negatives (in truth but not called)

        # For calculating TN, we need to estimate total possible sites
        # This is tricky for variant calling - using final_variants as proxy for called sites
        # and truth_total as proxy for true variant sites
        Total_called = final_variants,
        Total_truth = truth_total,

        # Standard metrics
        Sensitivity = TP / (TP + FN), # Recall/TPR
        Precision = TP / (TP + FP), # Precision/PPV
        F1_score = 2 * (Precision * Sensitivity) / (Precision + Sensitivity),

        # False Discovery Rate
        FDR = FP / (TP + FP) # False Discovery Rate

        # # False Positive Rate (requires estimate of total negatives)
        # # In variant calling context, this is often approximated as:
        # FPR = FP / (FP + TP), # Simplified FPR for variant calling

        # # Balanced accuracy (requires specificity)
        # # For variant calling, we can use a simplified version:
        # balanced_accuracy = (sensitivity + precision) / 2, # Simplified balanced accuracy

        # # Specificity is challenging in variant calling without knowing total genome sites
        # # Using a proxy based on called variants
        # specificity = ifelse(total_called - TP > 0,
        #     (total_called - TP - FP) / (total_called - TP),
        #     NA
        # )
    )  |> 
    # Handle NaN and Inf values
    mutate_if(is.numeric, ~ ifelse(is.nan(.) | is.infinite(.), 0, .)) |> 
    # Round to 4 decimal places
    mutate_if(is.numeric, ~ round(., 4))

# Create the final performance table with requested columns
performance_table <- performance_metrics |>
    rename(Mutect2_Call_Depth = depth) |>
    select(
        Mutect2_Call_Depth, TP, FP, FN,
        Total_called, Total_truth, Sensitivity, Precision, F1_score, FDR
    ) |> 
    mutate(
        Sensitivity = paste0(
            round(Sensitivity * 100, 2), "%"
        ),
        Precision = paste0(
            round(Precision * 100, 2), "%"
        ),
        F1_score = paste0(
            round(F1_score * 100, 2), "%"
        ),
        FDR = paste0(
            round(FDR * 100, 2), "%"
        )
    )

# Create a nicely formatted table
print(performance_table)

# Save the table to a CSV file
write_csv(
    performance_table, 
    file.path(work_dir, "mutect2_performance_metrics.csv")
)

# Save the table as plot
tg = gridExtra::tableGrob(performance_table)
h = grid::convertHeight(sum(tg$heights), "in", TRUE)
w = grid::convertWidth(sum(tg$widths), "in", TRUE)

ggsave(
    here(work_dir, "mutect2_performance_metrics.pdf"), 
    tg, 
    width=w, 
    height=h
)

# formatted_table <- knitr::kable(performance_table,
#     format = "html",
#     caption = "Mutect2 Performance Metrics Across Different Depths",
#     digits = 4
# ) %>%
#     kableExtra::kable_styling(
#         bootstrap_options = c("striped", "hover", "condensed"),
#         full_width = FALSE
#     )
# print(formatted_table)
# hcc1395 <- performance_table
# na12878 <- performance_table

compare_table <- bind_rows(
    na12878,
    hcc1395
) |> 
    mutate(
        truth_dataset = c("NA12878", "HCC1395")
    ) |> 
    relocate(
        truth_dataset, .before = Mutect2_Call_Depth
    )

tg = gridExtra::tableGrob(compare_table)
h = grid::convertHeight(sum(tg$heights), "in", TRUE)
w = grid::convertWidth(sum(tg$widths), "in", TRUE)

ggsave(
    here(work_dir, "benchmark_performance_metrics.pdf"), 
    tg, 
    width=w, 
    height=h
)
