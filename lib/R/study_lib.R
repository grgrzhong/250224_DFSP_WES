SavePlot <- function(
    plot, width = 8, height = 6, only_png = TRUE, dir, filename
) {

    ### Save png or pdf plots using ggplot2

    fs::dir_create(here(dir))

    if (only_png) {
        ### Save only png
        img_type <- ".png"
        ggsave(
            filename = here(dir, paste0(filename, img_type)),
            plot = plot,
            width = width,
            height = height,
            units = "in",
            dpi = 300,
            device = ifelse(img_type == ".png", png, cairo_pdf)
        )

    } else {
        ### Save both png and pdf
        for (img_type in c(".png", ".pdf")) {
            ggsave(
                filename = here(dir, paste0(filename, img_type)),
                plot = plot,
                width = width,
                height = height,
                units = "in",
                dpi = 300,
                device = ifelse(img_type == ".png", png, cairo_pdf)
            )
        }
    }
}

figure_theme <- function(line_width = 0.3, base_size = 8) {
    ## General theme for publication figures
    # theme_classic(base_size = 10, base_family = "Arial") +
    theme_classic(
        base_size = base_size,
        base_family = "Arial",
        base_line_size = line_width
    ) +
        theme(
            plot.title = element_text(
                face = "plain", #"bold", 
                hjust = 0.5, size = base_size
            ),
            strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "right",
            legend.title = element_text(
                hjust = 0,
                size = base_size - 1
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = base_size -1),
            legend.box = "vertical",
            # legend.box = "horizontal",
            legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                linewidth = line_width, color = "black"
            ),
            # axis.line = element_line(color = "black", linewidth = line_width),
            # axis.line = element_blank(),
            # panel.border = element_rect(
            #     linewidth = line_width, color = "black", fill = NA
            # ),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        )
}

figure_theme2 <- function() {
    ## General theme for publication figures
    # theme_classic(base_size = 10, base_family = "Arial") +
    line_width <- 0.3
    base_size <- 8
    theme_bw(
        base_size = base_size,
        base_family = "Arial",
        base_line_size = line_width
    ) +
        theme(
            plot.title = element_text(face = "plain", hjust = 0.5, size = 8),
            # strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "right",
            legend.title = element_text(
                hjust = 0,
                size = 8
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = base_size -1),
            legend.box = "vertical",
            # legend.box = "horizontal",
            legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                linewidth = line_width, #color = "black"
            ),
            # axis.line = element_line(color = "black", linewidth = line_width),
            # axis.line = element_blank(),
            # panel.border = element_rect(
            #     linewidth = line_width, color = "black", fill = NA
            # ),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        )
}

plot_theme <- function() {
    # theme_classic(base_size = 10, base_family = "Arial") +
    # theme_classic() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
            # strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "top",
            legend.title = element_text(
                hjust = 0, 
                size = 8
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 8),
            legend.box = "vertical",
            # legend.box = "horizontal",
            legend.box.just = "left", # Align legend box to the left,
            text = element_text(size = 11),
            axis.text = element_text(color = "black"),
            # axis.ticks = element_line(linewidth = 0.3, color = "black"),    # Set the thickness of axis ticks
            # axis.line = element_line(color = "black", linewidth = 0.5),
            # axis.line = element_blank(),
            # axis.line = element_line(linewidth = 0.5, color = "black"),     # Set the thickness of axis lines
            # panel.border = element_rect(linewidth = 0.5, color = "black", fill = NA),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), units = "cm")
        )

}

TestVariantFilter <- function(
    filter_params,
    maf,
    is_save = FALSE,
    fig_dir = "figures/annovar",
    fig_name = "compare_filter"
) {
    
    # Apply filtering with parameters
    maf_filtered <- maf |>
        # Filter by sequencing quality
        filter(
            DP >= filter_params$min_depth | is.na(DP),
            (is.na(AF) | as.numeric(AF) >= filter_params$min_vaf)
        ) |>
        # Filter by variant impact
        filter(
            !Variant_Classification %in% filter_params$exclude_classifications,
            !ExonicFunc.refGene %in% filter_params$exclude_exonic_func
        ) |>
        # Filter population variants
        # # Filter by pathogenicity predictions
        # filter(
        #   SIFT_pred == "D" |
        #   Polyphen2_HDIV_pred %in% c("D", "P") |
        #   MutationTaster_pred %in% c("D", "A") |
        #   (!is.na(CADD_phred) & as.numeric(CADD_phred) >= filter_params$min_cadd_phred) |
        #   is.na(SIFT_pred)
        # )
        mutate(
            gnomAD_exome_ALL_num = suppressWarnings(
                as.numeric(gnomAD_exome_ALL)
            )
        ) |>
        filter(
            is.na(gnomAD_exome_ALL_num) |
                gnomAD_exome_ALL_num <= filter_params$max_population_freq
        ) |>
        select(-gnomAD_exome_ALL_num)

    # Count variants per sample after filtering
    variants_per_sample_after <- maf_filtered |>
        count(Tumor_Sample_Barcode, name = "n_variants_after")

    # Merge before and after counts
    variant_comparison <- variants_per_sample_before |>
        left_join(variants_per_sample_after, by = "Tumor_Sample_Barcode") |>
        mutate(
            n_variants_after = replace_na(n_variants_after, 0),
            percent_retained = round(n_variants_after / n_variants_before * 100, 1)
        )

    # Print summary statistics
    cat(
        "Total variants after filtering:",
        sum(variant_comparison$n_variants_after), "\n"
    )

    cat(
        "Average variants per sample before:",
        mean(variant_comparison$n_variants_before), "\n"
    )

    cat(
        "Average variants per sample after:",
        mean(variant_comparison$n_variants_after), "\n"
    )

    cat(
        "Average percentage retained:",
        mean(variant_comparison$percent_retained), "%\n"
    )

    # Create comparison visualizations
    ## 1. Bar plot of variant counts before/after by sample
    p <- variant_comparison |>
        pivot_longer(
            cols = c(n_variants_before, n_variants_after),
            names_to = "filter_status",
            values_to = "variant_count"
        ) |>
        mutate(
            filter_status = ifelse(
                filter_status == "n_variants_before",
                "Before filtering", "After filtering"
            ),
            filter_status = factor(filter_status, levels = c("Before filtering", "After filtering"))
        ) |>
        ggplot(aes(x = variant_count, fill = filter_status)) +
        geom_histogram(
            binwidth = 500, alpha = 0.7, position = "identity", boundary = 0
        ) +
        facet_wrap(~filter_status, ncol = 1, scales = "free_y") +
        scale_fill_manual(
            values = c("Before filtering" = "#1f77b4", "After filtering" = "#ff7f0e")
        ) +
        # Add mean lines with labels
        geom_vline(
            data = . %>% group_by(filter_status) %>%
                summarize(mean_count = mean(variant_count)),
            aes(xintercept = mean_count),
            linetype = "dashed", size = 1, color = "red"
        ) +
        # Add mean value labels
        geom_text(
            data = . %>% group_by(filter_status) %>%
                summarize(mean_count = mean(variant_count)),
            aes(x = mean_count, y = Inf, label = sprintf("Mean: %d", round(mean_count))),
            hjust = -0.1, vjust = 1.5, color = "black", size = 7, size.unit = "pt"
        ) +
        # Colors and theme
        labs(
            title = "Distribution of variant counts per sample",
            x = "Number of variants",
            y = "Number of samples",
            fill = "Filter status"
        ) +
        figure_theme() +
        theme(
            # strip.text = element_text(size = 12, face = "bold"),
            # axis.title = element_text(size = 11),
            legend.position = "none"
        )

    if (is_save) {
    
        SavePlot(
            plot = p,
            width = 4,
            height = 4,
            dir = fig_dir,
            filename = fig_name
        )
    }

    list(
            maf_filtered = maf_filtered,
            variant_counts = variant_comparison,
            plot = p
        )
}