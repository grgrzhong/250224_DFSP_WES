# Load required libraries
library(VariantAnnotation)
library(ggbeeswarm)
source(here::here("bin/R/lib/study_lib.R"))

# Mutect2 with balcklist and RepeatMasker filtering, tumour normal
vcf_dir1 <- here("data/wes/mutect2")

# Mutect2 with balcklist and RepeatMasker filtering, tumour only
# vcf_dir2 <- here("data/wes/mutect2_tumour_only")

# Mutect2 without blacklisting and RepeatMasker filtering, tumour normal
vcf_dir2 <- here("data/wes/mutect2_old")

# group_name <- c("Mutect2_Tumour_Normal", "Mutect2_Tumour_Only")
group_name <- c("With_Black_RepeatMasker_filter", "Without_Black_RepeatMasker_filter")

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
vcf_files1 <- dir_ls(vcf_dir1, recurse = TRUE, glob = "*final.vcf.gz")
vcf_files2 <- dir_ls(vcf_dir2, recurse = TRUE, glob = "*.final.vcf.gz")

# Extract sample names
get_sample_name <- function(file_path) {
    file_path |>
        path_file() |>
        str_remove("\\.vcf(\\.gz)?$") |>
        str_remove(".final")
}

samples1 <- vcf_files1 |> map_chr(get_sample_name)
samples2 <- vcf_files2 |> map_chr(get_sample_name)

# Find common samples
common_samples <- intersect(samples1, samples2)
message("Found ", length(common_samples), " common samples\n")

# Create sample file mapping
sample_mapping <- tibble(
    sample = common_samples,
    file1 = map_chr(
        common_samples,
        \(s) vcf_files1[str_detect(vcf_files1, s)][1]
    ),
    file2 = map_chr(
        common_samples,
        \(s) vcf_files2[str_detect(vcf_files2, s)][1]
    )
)

# Process each sample with better error handling and clearer steps
results1 <- sample_mapping |>
    dplyr::select(sample, file1) |>
    mutate(
        count = map2(
            file1,
            sample,
            \(f, s) {
                cat("Processing:", s, "\n")
                count_variants(f)
            },
            .progress = TRUE
        )
    ) |>
    unnest(count) |>
    dplyr::rename(
        data1_total = total_variants,
        data1_pass = pass_variants
    ) |> 
    ## We keep the total as it is same as the pass variants
    dplyr::select(-file1)

# Step 2: Process filtered files -----
results2 <- sample_mapping |>
    dplyr::select(sample, file2) |>
    mutate(
        count = map2(
            file2,
            sample,
            \(f, s) {
                cat("Processing:", s, "\n")
                count_variants(f)
            },
            .progress = TRUE
        )
    ) |>
    unnest(count) |>
    dplyr::rename(
        data2_total = total_variants,
        data2_pass = pass_variants
    ) |> 
    ## We keep the total as it is same as the pass variants
    dplyr::select(-file2)

# Step 3: Combine results ----
results <- results1 |> 
    left_join(results2, by = "sample") |>
    dplyr::select(-matches("pass")) |>
    ## Remove pass columns as they are not needed for this analysis
    mutate(
        difference = data1_total - data2_total
    )

# Distribution of total variants Before filtering ------------

# Function to plot the distribution of total variants
compare_variants <- function(
    results,
    group_name,
    title
) {

    ## Get the mean and median data
    # mean_data1 <- round(mean(results$data1_total, na.rm = TRUE), 0)
    # median_data1 <- round(median(results$data1_total, na.rm = TRUE),0)
    # mean_data2 <- round(mean(results$data2_total, na.rm = TRUE), 0)
    # median_data2 <- round(median(results$data2_total, na.rm = TRUE),0)
    
    ## Prepare the data for plotting
    plot_data <- results |> 
        select(-difference) |>
        pivot_longer(
            cols = c(data1_total, data2_total),
            names_to = "group",
            values_to = "total_variants"
        ) |> 
        mutate(
            group = case_when(
                group == "data1_total" ~ group_name[1],
                group == "data2_total" ~ group_name[2]
            ) |> 
            factor(
                levels = group_name,
                ordered = TRUE
            )
        )

    ## Calculate mean and median for each group
    summary_stats <- plot_data |>
        group_by(group) |>
        summarise(
            mean_val = mean(total_variants, na.rm = TRUE),
            median_val = median(total_variants, na.rm = TRUE),
            .groups = "drop"
        )

        # annotate(
        #     geom = "text",
        #     x = I(0.5),
        #     y = I(0.8),
        #     label = paste0("Mean: ", mean_data1),
        #     hjust = 0,
        #     color = "red", 
        #     size.unit = "pt",
        #     size = 8
        # ) +
    text_size <- 8

    plot <- plot_data |>
        ggplot(aes(x = total_variants, fill = group)) +
        geom_histogram(bins = 30, alpha = 0.7, color = "white") +
        facet_wrap(~ group, ncol = 1) +
        geom_vline(
            data = summary_stats,
            aes(xintercept = mean_val),
            color = "red", linetype = "dashed", linewidth = 0.5
        ) +
        geom_vline(
            data = summary_stats,
            aes(xintercept = median_val),
            color = "darkred", linetype = "dotted", linewidth = 0.5
        ) +
        geom_text(
            data = summary_stats,
            aes(x = Inf, y = Inf, label = paste0("Mean: ", round(mean_val, 0))),
            hjust = 1.1, vjust = 2,
            color = "red", size = text_size, size.unit = "pt"
        ) +
        geom_text(
            data = summary_stats,
            aes(x = Inf, y = Inf, label = paste0("Median: ", round(median_val, 0))),
            hjust = 1.1, vjust = 3.5,
            color = "darkred", size = text_size, size.unit = "pt"
        ) +
        labs(
            title = title,
            x = "Number of Variants",
            y = "Number of Samples"
        ) +
        theme_minimal() +
        plot_theme() +
        theme(legend.position = "none")
    
    plot
}

# Set output directory
# out_dir <- here("results/mutect2_normal")
# filename <- "mutect2_normal_vs_nonormal"
# title = "Distribution of total variants\nMutect2 with vs without Normal"

out_dir <- here("results/black_repeatmasker")
filename <- "black_repeatmasker_with_vs_without"
title = "Distribution of total variants\nBlack list & repeatmasker with vs without filtering"

SavePlot(
    plot = compare_variants(results, group_name, title),
    width = 3.5, 
    height = 4, 
    dir = out_dir, 
    filename = filename
)

## Save the results to a CSV file
export <- results |> 
    rename(
        with_normal_total = data1_total,
        without_normal_total = data2_total,
        total_difference = difference
    )

write_excel_csv(
    export,
    file = here(out_dir, paste0(filename, ".csv"))
)

## Blacklist and RepeatMasker filtering ------------
out_dir <- here("results/black_repeatmasker")
filename <- "black_repeatmasker_with_vs_without"
title = "Distribution of total variants\nBlack list & repeatmasker with vs without filtering"

SavePlot(
    plot = compare_variants(results, group_name, title),
    width = 3.5, 
    height = 4, 
    dir = out_dir, 
    filename = filename
)

## Save the results to a CSV file
export <- results |> 
    rename(
        with_black_repeatmasker_filter = data1_total,
        without_black_repeatmasker_filter = data2_total,
        total_difference = difference
    )

write_excel_csv(
    export,
    file = here(out_dir, paste0(filename, ".csv"))
)
