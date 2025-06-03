## Load required libraries
suppressPackageStartupMessages(
    suppressWarnings({
        library(maftools)
        library(qs)
        library(fs)
        library(here)
        library(readxl)
        library(writexl)
        library(tidyverse)
    })
)

SavePlot <- function(
    plot, width = 8, height = 6, only_png = FALSE, dir, filename) {
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
        
        message(
            "Saved plot: ", here(dir, paste0(filename, img_type))
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
            
            message(
                "Saved plot: ", here(dir, paste0(filename, img_type))
            )
        }
    }
}

plot_theme <- function(line_width = 0.3, base_size = 8) {
    ## General theme for publication figures
    # theme_classic(base_size = 10, base_family = "Arial") +
    # theme_classic(
    #     base_size = base_size,
    #     base_family = "Arial",
    #     base_line_size = line_width
    # ) +
    theme(
            plot.title = element_text(
                face = "plain", # "bold",
                # hjust = 0.5, 
                size = base_size
            ),
            # strip.background = element_blank(),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            # legend.position = "right",
            legend.title = element_text(
                # hjust = 0,
                size = base_size - 1
                # face = "bold"
            ),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = base_size - 1),
            # legend.box = "vertical",
            # legend.box = "horizontal",
            # legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            # axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                # color = "black",
                linewidth = line_width
            ),
            # axis.line = element_line(color = "black", linewidth = line_width),
            # axis.line = element_blank(),
            # panel.border = element_rect(
            #     linewidth = line_width, color = "black", fill = NA
            # ),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), units = "cm")
        )
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
                face = "plain", # "bold",
                # hjust = 0.5, 
                size = base_size
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
            legend.text = element_text(size = base_size - 1),
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
            legend.text = element_text(size = base_size - 1),
            legend.box = "vertical",
            # legend.box = "horizontal",
            legend.box.just = "left", # Align legend box to the left,
            # text = element_text(size = 11),
            axis.text = element_text(color = "black"),
            # panel.grid = element_blank(),
            axis.ticks = element_line(
                linewidth = line_width, # color = "black"
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
    maf) {
    # Apply filtering with parameters
    maf_filtered <- maf |>
        ## Filter by sequencing quality
        filter(
            DP >= filter_params$min_depth | is.na(DP),
            (is.na(AF) | as.numeric(AF) >= filter_params$min_vaf),
            (is.na(VAD) | as.numeric(VAD) >= filter_params$min_vad),
        ) |>
        ## Filter by variant impact
        filter(
            !Variant_Classification %in% filter_params$exclude_classifications,
            !ExonicFunc.refGene %in% filter_params$exclude_exonic_func,
            (Func.refGene %in% filter_params$include_func | is.na(Func.refGene))
        ) |>
        ## Filter population variants
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
    post_filter_count <- maf_filtered |>
        count(Tumor_Sample_Barcode, name = "n_variants_after")

    # Merge before and after counts
    filter_res <- post_filter_count |>
        mutate(
            n_variants_after = replace_na(n_variants_after, 0),
            # percent_retained = round(n_variants_after / n_variants_before * 100, 1)
        )

    list(
        maf_filtered = maf_filtered,
        filter_res = filter_res
    )
}

PlotVariantFilter <- function(
    is_save = FALSE,
    fig_dir = "figures/variant_filtering",
    fig_name = "filter_params") {
    # Create comparison visualizations
    ## 1. Bar plot of variant counts before/after by sample
    p <- filter_res |>
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
            binwidth = 200, alpha = 0.7, position = "identity", boundary = 0
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
            linetype = "dashed", linewidth = 0.5, color = "red"
        ) +
        # Add mean value labels
        geom_text(
            data = . %>% group_by(filter_status) %>%
                summarize(mean_count = mean(variant_count)),
            aes(x = mean_count, y = Inf, label = sprintf("Mean: %d", round(mean_count))),
            hjust = -0.1, vjust = 1.5, color = "black", size = 7, size.unit = "pt"
        ) +
        # Add filter parameters as annotation only on "After filtering" panel
        geom_text(
            data = data.frame(
                filter_status = "Before filtering",
                x = 20000, # Adjust this position as needed
                y = I(0.3) # Use relative position on y-axis
            ),
            aes(x = x, y = y, label = sprintf(
                "Filter params:\nread depth >= %d\nvaf >= %.2f\nvad >= %d\npop_freq <= %.4f",
                filter_params$min_depth,
                filter_params$min_vaf,
                filter_params$min_vad,
                filter_params$max_population_freq
            )),
            hjust = 0, # Left-aligned text
            vjust = 0, # Top-aligned text
            size = 6,
            size.unit = "pt",
            color = "black",
            fontface = "plain"
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
            height = 3,
            dir = fig_dir,
            filename = fig_name
        )
    }
}

MafOncoPlot <- function(
    maf,
    genes = NULL,
    top_n_genes = 30,
    clinicalFeatures = NULL,
    annotationColor = NULL,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = TRUE,
    removeNonMutated = FALSE,
    titleText = "OncoPlot",
    fontSize = 0.7,
    width = 10,
    height = 8,
    fig_dir = "figures/oncoplot",
    fig_name = "oncoplot") {
    # Create directory if it doesn't exist
    fs::dir_create(here(fig_dir))

    # If no genes specified, get top mutated ones
    if (is.null(genes)) {
        gene_summary <- getGeneSummary(maf)
        n_genes <- min(top_n_genes, nrow(gene_summary))

        # Check if we have enough genes
        if (n_genes < 2) {
            warning("Not enough genes for oncoplot (minimum 2 required)")
            return(NULL)
        }

        genes <- gene_summary[1:n_genes, "Hugo_Symbol"]
        genes <- genes$Hugo_Symbol
    }

    # Save PDF version
    pdf_file <- here(fig_dir, paste0(fig_name, ".pdf"))
    pdf(pdf_file, width = width, height = height)

    p <- oncoplot(
        maf = maf,
        genes = genes,
        clinicalFeatures = clinicalFeatures,
        annotationColor = annotationColor,
        sortByAnnotation = sortByAnnotation,
        showTumorSampleBarcodes = showTumorSampleBarcodes,
        removeNonMutated = removeNonMutated,
        fontSize = fontSize,
        titleText = titleText
    )

    dev.off()
    message(paste("Saved plots to:", pdf_file))

    # save PNG version
    png_file <- here(fig_dir, paste0(fig_name, ".png"))
    png(png_file, width = width, height = height, res = 300, units = "in")

    oncoplot(
        maf = maf,
        genes = genes,
        clinicalFeatures = clinicalFeatures,
        annotationColor = annotationColor,
        sortByAnnotation = sortByAnnotation,
        showTumorSampleBarcodes = showTumorSampleBarcodes,
        removeNonMutated = removeNonMutated,
        fontSize = fontSize,
        titleText = titleText
    )
    dev.off()

    message(paste("Saved plots to:", png_file))

    print(p)
}

LoadSampleInfo <- function() {
    read_xlsx(
        here(
            "data/clinical/DFSP_multiomics_sample_list_updated_20250509.xlsx"
        )
    )
}

AddCancerHotspot <- function(
    maf,                           # maf data in dataframe format
    hotspot = NULL,                # path to the hotspot file
    qvalue = NULL,                 # qvalue threshold
    median_allele_freq_rank = NULL, # median Allele Frequency Rank threshold 
    log10_pvalue = NULL           # log10 pvalue threshold

) {
    if (is.null(hotspot)) {

        snv_hotspots <- read_xlsx(
            here("data/clinical/hotspots_v2.xlsx"),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here("data/clinical/hotspots_v2.xlsx"),
            sheet = "INDEL-hotspots"
        ) 
        
    } else {
        
        snv_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "INDEL-hotspots"
        )
    }

    snv_hotspots <- snv_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            ref_aa = str_extract(Reference_Amino_Acid, "^[A-Z*]"),
            pos_aa = as.character(Amino_Acid_Position),
            var_aa = str_extract(Variant_Amino_Acid, "^[A-Z*]"),
            aaChange = paste0("p.", ref_aa, pos_aa, var_aa)
        ) |>
        distinct() |>
        mutate(
            snv_hotspot = paste(
                Hugo_Symbol, aaChange,
                sep = "_"
            )
        )
    
    indel_hotspots <- indel_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            aaChange = str_extract(Variant_Amino_Acid, "^[^:]+")
        ) |>
        mutate(
            aaChange = paste0("p.", aaChange)
        ) |>
        distinct() |>
        mutate(
            indel_hotspot = paste0(
                Hugo_Symbol, "_", aaChange
            )
        )
    
    ## Match the variant info from the aaChange column in maf and
    ## Variant_Amino_Acid column in the hotspot data
    hotspots <- list(
        snv_hotspots = snv_hotspots,
        indel_hotspots = indel_hotspots
    )

    message(
        sprintf(
            "Total SNV hotspots: %d", nrow(hotspots$snv_hotspots)
        ), "\n",
        sprintf(
            "Total INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
        )
    )

    ## Apply filtering only if filter parameters are provided
    if (
        !is.null(qvalue) || !is.null(median_allele_freq_rank) ||
        !is.null(log10_pvalue)
    ) {
        ## build filter conditions
        filter_expr <- list()
    
        if (!is.null(qvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(qvalue < !!qvalue)
            )
        }
        
        if (!is.null(median_allele_freq_rank)) {
            filter_expr <- append(
                filter_expr, 
                rlang::expr(Median_Allele_Freq_Rank > !!median_allele_freq_rank)
            )
        }
        
        if (!is.null(log10_pvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(log10_pvalue > !!log10_pvalue)
            )
        }

        ## Apply filtering to select the high confidence hotspot variants
        hotspots <- map(
            hotspots,
            ~ .x |> filter(!!!filter_expr)
        )

        message("Applied filters:")
        if (!is.null(qvalue)) message("  - qvalue < ", qvalue)
        if (!is.null(median_allele_freq_rank)) message("  - Median_Allele_Freq_Rank > ", median_allele_freq_rank)
        if (!is.null(log10_pvalue)) message("  - log10_pvalue > ", log10_pvalue)

        message(
            sprintf(
                "After filtering SNV hotspots: %d", nrow(hotspots$snv_hotspots)
            ), "\n",
            sprintf(
                "After filtering INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
            )
        )

    } else {

        message("No filters applied - using all hotspots")
    }

    all_hotspots <- c(
        hotspots$snv_hotspots |> pull(snv_hotspot),
        hotspots$indel_hotspots |> pull(indel_hotspot)
    )

    maf_hotspot <- maf |>
        # select(
        #     Tumor_Sample_Barcode, Variant_Classification, Hugo_Symbol, aaChange
        # ) |>
        # filter(
        #     # Variant_Classification == "Missense_Mutation"
        #     Variant_Classification == "Nonsense_Mutation"
        # ) |>
        mutate(match_change = paste0(Hugo_Symbol, "_", aaChange)) |>
        mutate(
            is_hotspot = if_else(
                match_change %in% all_hotspots,
                TRUE,
                FALSE
            )
        ) |>
        arrange(desc(is_hotspot)) |>
        select(-match_change) |>
        relocate(
            is_hotspot,
            .after = aaChange
        )
    
    is_hotspot_data <- maf_hotspot |>
        filter(is_hotspot) |>
        select(
            Tumor_Sample_Barcode, Variant_Classification, Hugo_Symbol, aaChange,
            is_hotspot
        )
    
    message(
        sprintf(
            "Number of hotspot variants in the maf: %d", nrow(is_hotspot_data)
        )
    )

    return(maf_hotspot)
}

LoadCancerHotspot <- function(
    hotspot = NULL,                # path to the hotspot file
    qvalue = NULL,                 # qvalue threshold
    median_allele_freq_rank = NULL, # median Allele Frequency Rank threshold 
    log10_pvalue = NULL           # log10 pvalue threshold

) {
    if (is.null(hotspot)) {

        snv_hotspots <- read_xlsx(
            here("data/clinical/hotspots_v2.xlsx"),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here("data/clinical/hotspots_v2.xlsx"),
            sheet = "INDEL-hotspots"
        ) 
        
    } else {
        
        snv_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "SNV-hotspots"
        )
        

        indel_hotspots <- read_xlsx(
            here(hotspot),
            sheet = "INDEL-hotspots"
        )
    }

    snv_hotspots <- snv_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            ref_aa = str_extract(Reference_Amino_Acid, "^[A-Z*]"),
            pos_aa = as.character(Amino_Acid_Position),
            var_aa = str_extract(Variant_Amino_Acid, "^[A-Z*]"),
            aaChange = paste0("p.", ref_aa, pos_aa, var_aa)
        ) |>
        distinct() |>
        mutate(
            snv_hotspot = paste(
                Hugo_Symbol, aaChange,
                sep = "_"
            )
        )
    
    indel_hotspots <- indel_hotspots |>
        select(
            Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, qvalue, Median_Allele_Freq_Rank
        ) |>
        mutate(
            aaChange = str_extract(Variant_Amino_Acid, "^[^:]+")
        ) |>
        mutate(
            aaChange = paste0("p.", aaChange)
        ) |>
        distinct() |>
        mutate(
            indel_hotspot = paste0(
                Hugo_Symbol, "_", aaChange
            )
        )
    
    ## Match the variant info from the aaChange column in maf and
    ## Variant_Amino_Acid column in the hotspot data
    hotspots <- list(
        snv_hotspots = snv_hotspots,
        indel_hotspots = indel_hotspots
    )

    message(
        sprintf(
            "Total SNV hotspots: %d", nrow(hotspots$snv_hotspots)
        ), "\n",
        sprintf(
            "Total INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
        )
    )

    ## Apply filtering only if filter parameters are provided
    if (
        !is.null(qvalue) || !is.null(median_allele_freq_rank) ||
        !is.null(log10_pvalue)
    ) {
        ## build filter conditions
        filter_expr <- list()
    
        if (!is.null(qvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(qvalue < !!qvalue)
            )
        }
        
        if (!is.null(median_allele_freq_rank)) {
            filter_expr <- append(
                filter_expr, 
                rlang::expr(Median_Allele_Freq_Rank > !!median_allele_freq_rank)
            )
        }
        
        if (!is.null(log10_pvalue)) {
            
            filter_expr <- append(
                filter_expr, rlang::expr(log10_pvalue > !!log10_pvalue)
            )
        }

        ## Apply filtering to select the high confidence hotspot variants
        hotspots <- map(
            hotspots,
            ~ .x |> filter(!!!filter_expr)
        )

        message("Applied filters:")
        if (!is.null(qvalue)) message("  - qvalue < ", qvalue)
        if (!is.null(median_allele_freq_rank)) message("  - Median_Allele_Freq_Rank > ", median_allele_freq_rank)
        if (!is.null(log10_pvalue)) message("  - log10_pvalue > ", log10_pvalue)

        message(
            sprintf(
                "After filtering SNV hotspots: %d", nrow(hotspots$snv_hotspots)
            ), "\n",
            sprintf(
                "After filtering INDEL hotspots: %d", nrow(hotspots$indel_hotspots)
            )
        )

    } else {

        message("No filters applied - using all hotspots")
    }

    hotspots
}


MergeAnnovarOutput <- function(
    annovar_dir,
    is_save = FALSE,
    save_dir = "data/wes/annotation/merged"
) {
    
    ## The annovar directory should contain the annovar output files
    annovar_dir <- here(annovar_dir)

    input_files <- dir_ls(annovar_dir, recurse = TRUE, glob = "*annovar.txt")
    
    ## Collect all annovar output files
    message(
        "Collecting ", length(input_files),
        " annovar output files from: ", annovar_dir
    )
    
    maf_data <- annovarToMaf(input_files) |>
        as_tibble() |>
        mutate(
            ## If the tumour sample has a matched normal sample, Otherinfo13 = normal allell data, Otherinfo14 = tumour allell data
            ## if the tumour sample does not have a matched normal sample, Otherinfo13 = tumour allell data, Otherinfo14 = NA
            is_paired_normal = if_else(
                is.na(Otherinfo14) | Otherinfo14 == "NA" | Otherinfo14 == "",
                FALSE,
                TRUE
            )
        )
    
    # maf_data |> 
    #     select(
    #         Tumor_Sample_Barcode, Hugo_Symbol, aaChange, 
    #         Otherinfo13, Otherinfo14, is_paired_normal,
    #         AD, AF, DP
    #     )

    ## Extract the allele data for paired samples
    maf_data_paired <- maf_data |> filter(is_paired_normal)
    maf_data_unpaired <- maf_data |> filter(!is_paired_normal)

    n_paired_samples <- length(unique(maf_data_paired$Tumor_Sample_Barcode))
    n_unpaired_samples <- length(unique(maf_data_unpaired$Tumor_Sample_Barcode))

    if (n_paired_samples > 0) {

        message(
            sprintf(
                "Extracting allele data for tumor-normal sample: %d", 
                n_paired_samples
            )
        )
        
        maf_data_paired <- maf_data_paired |> 
            ## Extract the tumour allele data
            mutate(
                tumor_AD = sapply(
                    str_split(Otherinfo14, ":"),
                    function(x) x[2]
                )
            ) |> 
            separate(
                tumor_AD, 
                into = c("tumor_RAD", "tumor_VAD"), 
                sep = ",", 
                remove = TRUE,

            ) |> 
            mutate(
                tumor_RAD = as.numeric(tumor_RAD),
                tumor_VAD = as.numeric(tumor_VAD)
            ) |> 
            mutate(
                tumor_AF = as.numeric(
                        sapply(
                        str_split(Otherinfo14, ":"), 
                        function(x) x[3]
                    )
                ),
                tumor_DP = as.numeric(
                        sapply(
                        str_split(Otherinfo14, ":"), 
                        function(x) x[4]
                    )
                )
            ) |> 
            ## Extract the normal allele data
            mutate(
                normal_AD = sapply(
                    str_split(Otherinfo13, ":"),
                    function(x) x[2]
                )
            ) |> 
            separate(
                normal_AD, 
                into = c("normal_RAD", "normal_VAD"), 
                sep = ",", 
                remove = TRUE,

            ) |> 
            mutate(
                normal_RAD = as.numeric(normal_RAD),
                normal_VAD = as.numeric(normal_VAD)
            ) |> 
            mutate(
                normal_AF = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[3]
                    )
                ),
                normal_DP = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[4]
                    )
                )
            )

    } else {

        maf_data_paired <- maf_data_paired |> 
            mutate(
                tumor_RAD = NA_real_,
                tumor_VAD = NA_real_,
                tumor_AF = NA_real_,
                tumor_DP = NA_real_,
                normal_RAD = NA_real_,
                normal_VAD = NA_real_,
                normal_AF = NA_real_,
                normal_DP = NA_real_
            )
    }

    if (n_unpaired_samples > 0) {

        message(
            sprintf(
                "Extracting allele data for tumor-only sample: %d", 
                n_unpaired_samples
            )
        )

        maf_data_unpaired <- maf_data_unpaired |> 
            mutate(
                tumor_AD = sapply(
                    str_split(Otherinfo13, ":"),
                    function(x) x[2]
                )
            ) |> 
            separate(
                tumor_AD, 
                into = c("tumor_RAD", "tumor_VAD"), 
                sep = ",", 
                remove = TRUE,

            ) |> 
            mutate(
                tumor_RAD = as.numeric(tumor_RAD),
                tumor_VAD = as.numeric(tumor_VAD)
            ) |>
            mutate(
                tumor_AF = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[3]
                    )
                ),
                tumor_DP = as.numeric(
                        sapply(
                        str_split(Otherinfo13, ":"), 
                        function(x) x[4]
                    )
                )
            ) |> 
            mutate(
                normal_RAD = NA_real_,
                normal_VAD = NA_real_,
                normal_AF = NA_real_,
                normal_DP = NA_real_
            )
    } else {
        
        maf_data_unpaired <- maf_data_unpaired |> 
            mutate(
                tumor_RAD = NA_real_,
                tumor_VAD = NA_real_,
                tumor_AF = NA_real_,
                tumor_DP = NA_real_,
                normal_RAD = NA_real_,
                normal_VAD = NA_real_,
                normal_AF = NA_real_,
                normal_DP = NA_real_
            )
    }

    ## Combine the paired and unpaired samples
    all_columns <- colnames(maf_data_paired)
    
    maf_data_unpaired <- maf_data_unpaired |> 
        select(all_of(all_columns))
    
    stopifnot(
        all.equal(
            colnames(maf_data_paired),
            colnames(maf_data_unpaired)
        )
    )

    maf_tbl <- bind_rows(
        maf_data_paired,
        maf_data_unpaired
    ) |>
        mutate(
            gnomAD_exome_ALL = as.numeric(
                replace(gnomAD_exome_ALL, gnomAD_exome_ALL == ".", NA)
            )
        ) |> 
        ## reange the columns
        relocate(
            is_paired_normal, 
            Otherinfo13, Otherinfo14,
            AD, AF, DP, 
            tumor_RAD, tumor_VAD, tumor_AF, tumor_DP,
            normal_RAD, normal_VAD, normal_AF, normal_DP,
            .after = last_col()
        )

    ## Save the data
    if (is_save && !is.null(save_dir)) {
        # Save the merged maf file
        file_name <- "annovar_maf_merged.qs"
        
        message(
            "Saving merged annovar data: ", here(save_dir, file_name)
        )

        fs::dir_create(here(save_dir))
        
        qsave(maf_tbl, here(save_dir, file_name))
    }
    
    # Return a tibble
    maf_tbl

}

LoadMergedAnnovar <- function(annovar_dir = "data/wes/annotation/merged") {

    # Load the merged maf file
    file_name <- "annovar_maf_merged.qs"

    merged_maf <- here(annovar_dir, file_name)
    
    if (!file.exists(merged_maf)) {
    
        stop("Merged maf file not found: ", merged_maf)
    }
    
    message("Loading merged maf file: ", merged_maf)

    # Load the merged maf file and Clean the AD column
    maf_tbl <- qread(merged_maf)

    # Return a tibble
    maf_tbl
}

FilterMergedMaf <- function(
    maf_tbl,
    # Filtering parameters by selected maf columns
    filter_params = list(
        # Minimum read depth
        DP = list(op = ">=", value = 8),
        # Minimum variant allele depth
        VAD = list(op = ">=", value = 4),
        # Minimum variant allele frequency
        AF = list(op = ">=", value = 0.05),
        # Maximum population frequency
        gnomAD_exome_ALL = list(op = "<=", value = 0.01),
        
        # Exclude these classifications
        Variant_Classification = list(op = "out", value = c("Silent")),
        
        # Exclude these functional annotations
        Func.refGene = list(op = "out", value = c("synonymous_SNV"))
    ),
    print_summary = TRUE
) {

    # Summary variants before filtering
    summary_before <- maf_tbl |>
        group_by(Tumor_Sample_Barcode) |>
        summarise(n_variants = n()) |>
        summarise(
            across(
                n_variants,
                list(
                    total = ~ sum(.x),
                    mean = ~ mean(.x, na.rm = TRUE),
                    median = ~ median(.x, na.rm = TRUE)
                ),
                .names = "{.fn}_variants_before"
            )
        )

    # # Glimplese the maf data
    # colnames(maf_tbl)
    # unique(maf_tbl$Variant_Classification)
    # unique(maf_tbl$ExonicFunc.refGene)
    # unique(maf_tbl$Func.refGene)

    # ## Functional prediction
    # unique(maf_tbl$SIFT_pred)
    # unique(maf_tbl$Polyphen2_HDIV_pred)
    # unique(maf_tbl$MutationTaster_pred)
    # unique(maf_tbl$MutationAssessor_pred)
    # unique(maf_tbl$FATHMM_pred)
    # unique(maf_tbl$PROVEAN_pred)
    # unique(maf_tbl$MetaSVM_pred)
    # unique(maf_tbl$`M-CAP_pred`)

    maf_filter <- maf_tbl

    # Dynamically apply filters based on the provided filter_params
    for (param in names(filter_params)) {
        
        filter_rule <- filter_params[[param]]
        
        op <- filter_rule$op
        
        value <- filter_rule$value
        
        # Apply numeric filters
        if (is.numeric(value)) {

            if (op == ">=") {

                maf_filter <- maf_filter |> 
                    dplyr::filter(
                        is.na(.data[[param]]) | .data[[param]] >= value
                    )

            } else if (op == "<=") {

                maf_filter <- maf_filter |> 
                    dplyr::filter(
                        is.na(.data[[param]]) | .data[[param]] <= value
                    )
            }
        }
        
        # Apply categorical filters
        if (is.character(value) || is.factor(value)) {
            
            if (op == "in") {
                
                maf_filter <- maf_filter |> 
                    dplyr::filter(.data[[param]] %in% value)

            } else if (op == "out") {

                maf_filter <- maf_filter |> 
                    dplyr::filter(!(.data[[param]] %in% value))
            }
        }
    }
    
    # maf_filtered <- maf_tbl |> 
    #     ## Filter out variants with low read depth, variant allele depth,
    #     dplyr::filter(
    #         is.na(DP) | DP >= min_rd,
    #         is.na(VAD) | as.numeric(VAD) >= min_vad,
    #         is.na(AF) | as.numeric(AF) >= min_vaf
    #     ) |> 
    #     ## Filter out variants with high population frequency
    #     dplyr::filter(
    #         is.na(gnomAD_exome_ALL) | gnomAD_exome_ALL <= max_pop_freq
    #     ) |> 
    #     arrange(gnomAD_exome_ALL) |> 
    #     ## Filter out silent (synonymous) variants
    #     dplyr::filter(
    #         (Variant_Classification %in% c("Silent")),
    #         (Func.refGene %in% c("synonymous_SNV")),
    #         (ExonicFunc.refGene %in% c("synonymous SNV"))

    #     ) |> 
    #     ## Inlcude Filter out non-exonic variants
    #     dplyr::filter(
    #         Func.refGene %in% c(
    #             "exonic", "splicing", "UTR5"
    #         )
    #     )

    summary_after <- maf_filter |>
        group_by(Tumor_Sample_Barcode) |>
        summarise(n_variants = n()) |>
        summarise(
            across(
                n_variants,
                list(
                    total = ~ sum(.x),
                    mean = ~ mean(.x, na.rm = TRUE),
                    median = ~ median(.x, na.rm = TRUE)
                ),
                .names = "{.fn}_variants_after"
            )
        )

    # Print summary if requested
    if (print_summary) {
        message("Variant filtering summary:")
        message("---------------------------------------------------")
        message(
            sprintf(
                "Total variants before filtering: %d",
                summary_before$total_variants_before
            )
        )
        message(
            sprintf(
                "Total variants after filtering:  %d",
                summary_after$total_variants_after
            )
        )
        message(
            sprintf(
                "\nMean variants before:            %.2f",
                summary_before$mean_variants_before
            )
        )
        message(
            sprintf(
                "Mean variants after:             %.2f",
                summary_after$mean_variants_after
            )
        )
        message(
            sprintf(
                "\nMedian variants before:          %.2f",
                summary_before$median_variants_before
            )
        )
        message(
            sprintf(
                "Median variants after:           %.2f",
                summary_after$median_variants_after
            )
        )
    }
    
    # Return the filtered maf
    return(maf_filter)
}
