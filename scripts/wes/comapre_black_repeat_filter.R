library(VariantAnnotation)
library(ggbeeswarm)
source(here::here("bin/R/lib/study_lib.R"))

# Set data directory
unfiltered_dir <- here("data/wes/variant_calling/mutect2_without_black_repeat_filter")

## mutect2 = mutect2_filter
filtered_dir <- here("data/wes/variant_calling/mutect2_with_black_repeat_filter_new")

# Set output directory
out_dir <- here("results/black_repeat_filter")
dir_create(out_dir)

# Function to count variants in VCF file
count_variants <- function(vcf_file) {
    if (!file.exists(vcf_file)) {
        return(
            tibble(total_variants = NA_integer_, pass_variants = NA_integer_)
        )
    }

    tryCatch(
        {
            # Suppress VCF header warnings for duplicate keys
            suppressWarnings({
                # Adjust genome reference as needed
                vcf <- readVcf(vcf_file, "hg38")
            })

            total_vars <- nrow(vcf)
            pass_vars <- filt(vcf) |>
                (\(x) x == "PASS")() |>
                sum(na.rm = TRUE)

            tibble(
                total_variants = total_vars,
                pass_variants = pass_vars
            )
        },
        error = function(e) {
            cat("Error processing file:", vcf_file, "-", e$message, "\n")
            tibble(total_variants = NA_integer_, pass_variants = NA_integer_)
        }
    )
}

# Get list of VCF files from both directories
unfiltered_files <- dir_ls(
    unfiltered_dir,
    recurse = TRUE,
    glob = "*normalized_filtered.vcf.gz"
)

filtered_files <- dir_ls(
    filtered_dir,
    recurse = TRUE,
    glob = "*.final.vcf.gz"
)

# Extract sample names
get_sample_name_unfilter <- function(file_path) {
    file_path |>
        path_file() |>
        str_remove("\\.vcf(\\.gz)?$") |>
        str_remove("_normalized_filtered")
}

get_sample_name_filter <- function(file_path) {
    file_path |>
        path_file() |>
        str_remove("\\.vcf(\\.gz)?$") |>
        str_remove(".final")
}

unfiltered_samples <- unfiltered_files |> 
    map_chr(get_sample_name_unfilter)

filtered_samples <- filtered_files |> 
    map_chr(get_sample_name_filter)

# Find common samples
common_samples <- intersect(unfiltered_samples, filtered_samples)
cat("Found", length(common_samples), "common samples\n")

# Create sample file mapping
sample_mapping <- tibble(
    sample = common_samples,
    unfiltered_file = map_chr(
        common_samples,
        \(s) unfiltered_files[str_detect(unfiltered_files, s)][1]
    ),
    filtered_file = map_chr(
        common_samples,
        \(s) filtered_files[str_detect(filtered_files, s)][1]
    )
)

# Process each sample with better error handling and clearer steps
# Step 1: Process unfiltered files ----
unfiltered_results <- sample_mapping |>
    dplyr::select(sample, unfiltered_file) |>
    mutate(
        unfiltered_counts = map(
            unfiltered_file,
            \(f) {
                cat("  Processing:", basename(f), "\n")
                count_variants(f)
            },
            .progress = "Unfiltered files"
        )
    ) |>
    unnest(unfiltered_counts) |>
    dplyr::rename(
        unfiltered_total = total_variants,
        unfiltered_pass = pass_variants
    ) |> 
    ## We keep the total as it is same as the pass variants
    dplyr::select(-unfiltered_file)

# Step 2: Process filtered files -----
cat("\nStep 2: Processing filtered files...\n")
filtered_results <- sample_mapping |>
    dplyr::select(sample, filtered_file) |>
    mutate(
        filtered_counts = map(
            filtered_file,
            \(f) {
                cat("  Processing:", basename(f), "\n")
                count_variants(f)
            },
            .progress = "Filtered files"
        )
    ) |>
    unnest(filtered_counts) |>
    dplyr::rename(
        filtered_total = total_variants,
        filtered_pass = pass_variants
    ) |>
    dplyr::select(-filtered_file)

# Step 3: Combine results ----
results <- unfiltered_results |>
    left_join(filtered_results, by = "sample") |>
    filter(!is.na(unfiltered_total), !is.na(filtered_total))

cat("Successfully processed", nrow(results), "samples\n")

# Calculate differences and percentages
results <- results |>
    ## Remove pass columns as they are not needed for this analysis
    dplyr::select(-matches("pass")) |>
    mutate(
        total_removed = unfiltered_total - filtered_total,
        total_removal_pct = if_else(unfiltered_total > 0,
            (total_removed / unfiltered_total) * 100,
            0
        )
    )

write_csv(
    results,
    file = here(out_dir, "black_repeat_filter_results.csv")
)

range(results$total_removal_pct)
mean(results$total_removal_pct)
median(results$total_removal_pct)

# Distribution of total variants Before filtering ------------
mean_unfiltered <- round(
    mean(results$unfiltered_total, na.rm = TRUE),
    0
)
median_unfiltered <- round(
    median(results$unfiltered_total, na.rm = TRUE),
    0
)

plot <- results |>
    ggplot(aes(x = unfiltered_total)) +
    geom_histogram(
        bins = 30, fill = "coral", alpha = 0.7, color = "white"
    ) +
    geom_vline(
        xintercept = mean_unfiltered,
        color = "red", linetype = "dashed", linewidth = 0.5
    ) +
    annotate(
        geom = "text",
        x = I(0.5),
        y = I(0.8),
        label = paste0("Mean: ", mean_unfiltered),
        hjust = 0,
        color = "red", 
        size.unit = "pt",
        size = 8
    ) +
    geom_vline(
        xintercept = median_unfiltered,
        color = "darkred", linetype = "dotted", linewidth = 0.5
    ) +
    annotate(
        geom = "text",
        x = I(0.5),
        y = I(0.6),
        label = paste0("Median: ", median_unfiltered),
        hjust = 0,
        color = "red", 
        size.unit = "pt",
        size = 8
    ) +
    labs(
        title = "Distribution of total variants\nwithout Blacklist and Repeat Masker Filtering",
        x = "Number of Variants",
        y = "Number of Samples"
    ) +
    theme_minimal()+
    figure_theme2()

SavePlot(
    plot = plot,
    width = 4, 
    height = 3, 
    dir = out_dir, 
    filename = "black_repeat_filter_unfiltered_distribution"
)

# Distribution of total variants After filtering ------------
mean_filtered <- round(
    mean(results$filtered_total, na.rm = TRUE),
    0
)

median_filtered <- round(
    median(results$filtered_total, na.rm = TRUE),
    0
)

plot <- results |>
    ggplot(aes(x = filtered_total)) +
    geom_histogram(
        bins = 30, fill = "coral", alpha = 0.7, color = "white"
    ) +
    geom_vline(
        xintercept = mean_filtered,
        color = "red", linetype = "dashed", linewidth = 0.5
    ) +
    annotate(
        geom = "text",
        x = I(0.5),
        y = I(0.8),
        label = paste0("Mean: ", mean_filtered),
        hjust = 0,
        color = "red", 
        size.unit = "pt",
        size = 8
    ) +
    geom_vline(
        xintercept = median_filtered,
        color = "darkred", linetype = "dotted", linewidth = 0.5
    ) +
    annotate(
        geom = "text",
        x = I(0.5),
        y = I(0.6),
        label = paste0("Median: ", median_filtered),
        hjust = 0,
        color = "red", 
        size.unit = "pt",
        size = 8
    ) +
    labs(
        title = "Distribution of total variants\nwith Blacklist and Repeat Masker Filtering",
        x = "Number of Variants",
        y = "Number of Samples"
    ) +
    theme_minimal()+
    figure_theme2()

SavePlot(
    plot = plot,
    width = 4, 
    height = 3, 
    dir = out_dir, 
    filename = "black_repeat_filter_filtered_distribution"
)


# # Generate summary statistics for barplot
# summary_stats <- results |>
#     summarise(
#         samples_processed = n(),
#         # Before filtering stats
#         unfiltered_mean = mean(unfiltered_total, na.rm = TRUE),
#         unfiltered_median = median(unfiltered_total, na.rm = TRUE),
#         # After filtering stats
#         filtered_mean = mean(filtered_total, na.rm = TRUE),
#         filtered_median = median(filtered_total, na.rm = TRUE),
#         # Removal stats
#         across(c(total_removed, total_removal_pct),
#             list(
#                 mean = \(x) mean(x, na.rm = TRUE),
#                 median = \(x) median(x, na.rm = TRUE),
#                 sd = \(x) sd(x, na.rm = TRUE)
#             ),
#             .names = "{.col}_{.fn}"
#         )
#     )

# # Create data for comparison plots
# comparison_data <- results |>
#     select(sample, unfiltered_total, filtered_total) |>
#     pivot_longer(
#         cols = c(unfiltered_total, filtered_total),
#         names_to = "filter_status",
#         values_to = "variant_count"
#     ) |>
#     mutate(
#         filter_status = case_when(
#             filter_status == "unfiltered_total" ~ "Before Filtering",
#             filter_status == "filtered_total" ~ "After Filtering"
#         ) |> 
#         factor(
#             levels = c("Before Filtering", "After Filtering"),
#             ordered = TRUE
#         )
#     )

# # Create visualizations ----
# # Plot 1: Boxplot with beeswarm dots comparison
# p1 <- comparison_data |>
#     ggplot(aes(x = filter_status, y = variant_count, color = filter_status)) +
#     geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#     geom_beeswarm(alpha = 0.6, size = 1.2) +
#     # Add mean points
#     stat_summary(
#         fun = mean, geom = "point", shape = 23, size = 3,
#         fill = "white", color = "black", stroke = 1
#     ) +
#     # Add mean and median text annotations
#     stat_summary(
#         fun = mean, geom = "text",
#         aes(label = paste("Mean:", round(after_stat(y)))),
#         vjust = -1.5, hjust = 0.5, size = 3.5, color = "black"
#     ) +
#     stat_summary(
#         fun = median, geom = "text",
#         aes(label = paste("Median:", round(after_stat(y)))),
#         vjust = -0.5, hjust = 0.5, size = 3.5, color = "black"
#     ) +
#     scale_color_manual(
#         values = c(
#             "Before Filtering" = "coral",
#             "After Filtering" = "steelblue"
#         )
#     ) +
#     labs(
#         title = paste0(
#             "Variant Count Comparison: Before vs After Filtering\n",
#             " (", nrow(results), " samples)"
#         ),
#         # subtitle = str_glue("Across {nrow(results)} samples"),
#         x = "Filter Status",
#         y = "Number of Variants",
#         color = "Filter Status"
#     ) +
#     theme_minimal() +
#     figure_theme2() +
#     theme(legend.position = "none")

# SavePlot(
#     plot = p1,
#     width = 4, 
#     height = 3, 
#     dir = out_dir, 
#     filename = "black_repeat_filter_comparison"
# )

