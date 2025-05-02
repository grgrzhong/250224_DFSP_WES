
source(here::here("lib/R/study_lib.R"))

##############################################################################
## Merge the annovar output ----------------------
##############################################################################
input_dir <- here("data/wes/annotation/annovar")

input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*annovar.txt")
# length(input_files)
maf <- annovarToMaf(input_files)

## Clean the sample names
maf <- maf |> 
    mutate(
        Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_annovar", "")
    )

qsave(
  maf, 
  here("data/wes/annotation/merged/annovar_maf_merged.qs")
)

##############################################################################
## Merge the funcotator output --------------
##############################################################################
# input_dir <- here("data/wes/annotation/funcotator")
# input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*annotated.tsv")

# input_data <- map2(
#   input_files, names(input_files),
#   ~ {
#     maf_filtered <- read_tsv(.x, show_col_types = FALSE)
#     maf_filtered$Tumor_Sample_Barcode <- basename(.y) |> str_remove("\\..*$")

#     return(maf_filtered)
#   }
# )

# maf <- merge_mafs(mafs = input_data, verbose = TRUE)

# qsave(maf, here("data/wes/annotation/merged/funcotator_maf_merged.qs"))

##############################################################################
## Filtering the variants --------------
##############################################################################
maf <- qread(here("data/wes/annotation/merged/annovar_maf_merged.qs")) |> 
    separate(AD, into = c("RAD", "VAD"), sep = ",", remove = TRUE) |> 
    as_tibble()

colnames(maf)   
unique(maf$Variant_Classification)
unique(maf$Func.refGene)
unique(maf$ExonicFunc.refGene)

# Define filtering thresholds to filter out potenital sequencing errors or 
# artifacts while retaining true variants
filter_params <- list(
    min_depth = 10,               # Minimum read depth
    min_vaf = 0.05,               # Minimum variant allele frequency
    min_vad = 5,                  # Minimum variant allele depth, 3-10
    max_population_freq = 0.001,  # Maximum population frequency
    # Variant classifications to exclude
    exclude_classifications = c(
    "Silent", "Intron"
    ),
    # Exonic functions to exclude
    exclude_exonic_func = c(
    "synonymous_SNV"
    ),
    # functions to include
    include_func = c(             
    "exonic", "splicing"
    )
)

# Create a grid of parameter combinations to test
params_grid <- expand.grid(
    min_vaf = c(0.01, 0.02, 0.03, 0.04, 0.05),
    min_vad = c(3:10),
    min_depth = c(8:15)
)

# Create a directory to store all filter test results
filter_res_dir <- here("results/variant_filtering")
dir_create(filter_res_dir)

filter_results <- list()

pre_filter_counts <- maf |>
    count(Tumor_Sample_Barcode, name = "n_variants_before")

mean_pre_filter <- mean(
    pre_filter_counts$n_variants_before, na.rm = TRUE
)
cat("Pre-filter mean variants: ", mean_pre_filter, "\n")

median_pre_filter <- median(
    pre_filter_counts$n_variants_before, na.rm = TRUE
)
cat("Pre-filter median variants: ", median_pre_filter, "\n")

# Loop through all parameter combinations
for (i in 1:nrow(params_grid)) {

    # Create a new parameter set based on the base parameters
    current_params <- filter_params
    
    # Update with the current grid values
    current_params$min_vaf <- params_grid$min_vaf[i]
    current_params$min_vad <- params_grid$min_vad[i]
    current_params$min_depth <- params_grid$min_depth[i]
    
    # Create a descriptive name for this parameter combination
    param_name <- paste0(
        "vaf", current_params$min_vaf, 
        "_vad", current_params$min_vad, 
        "_depth", current_params$min_depth
    )
    
    message(paste0("Testing parameter combination: ", param_name))
    
    # Run the filter test
    filter_res <- TestVariantFilter(
        filter_params = current_params,
        maf = maf
    )
    
    mean_post_filter <- mean(
        filter_res$filter_res$n_variants_after, na.rm = TRUE
    )

    cat("Post-filter mean variants: ", mean_post_filter, "\n")

    median_post_filter <- median(
        filter_res$filter_res$n_variants_after, na.rm = TRUE
    )

    cat("Post-filter median variants: ", median_post_filter, "\n")

    # Get summary metrics about this filter combination
    filter_summary <- tibble(
        param_name = param_name,
        min_vaf = current_params$min_vaf,
        min_vad = current_params$min_vad,
        min_depth = current_params$min_depth,
        total_variants = nrow(maf),
        total_variants_retain = nrow(filter_res$maf_filtered),
        total_percent_return = round(
            total_variants_retain / total_variants * 100, 2
        ),
        mean_pre_filter = mean_pre_filter,
        mean_post_filter = mean_post_filter,
        mean_percent_retain = round(
            mean_post_filter / mean_pre_filter * 100, 2
        ),
        median_pre_filter = median_pre_filter,
        median_post_filter = median_post_filter,
        median_percent_retain = round(
            median_post_filter / median_pre_filter * 100, 2
        )
    )
    
    # Append to a summary file
    if (i == 1) {
        write_csv(
            filter_summary, 
            file.path(filter_res_dir, "variant_filter_summary.csv")
        )

    } else {
        write_csv(
            filter_summary, 
            file.path(filter_res_dir, "variant_filter_summary.csv"), 
            append = TRUE
        )
    }

    # Save the filtering results
    filter_results[[param_name]] <- filter_res
}

## Find the optimal filtering parameters, target mean variants: 1000
## we should use the median as the target value, as it is less affected by outliers
optimal_res <- read_csv(
    here(filter_res_dir, "variant_filter_summary.csv"),
    show_col_types = FALSE
) |> 
    arrange(desc(median_post_filter)) |>
    filter(param_name == "vaf0.01_vad3_depth10") |>
    filter(
        mean_post_filter > 1000 & 
        median_post_filter > 500
    )

## vaf = 0.01, vad =3, depth = 10
maf_list <- list(
    maf = maf,
    maf_filtered = filter_results[["vaf0.01_vad3_depth10"]][["maf_filtered"]]
)

qsave(maf_list, here(filter_res_dir, "variant_filter_results.qs"))

# Create a MAF object from filtered data for further analysis
maf_filtered <- read.maf(maf = filter_res$maf_filtered)

# Save filtered data
# qsave(
#   maf_filtered, 
#   here("data/wes/annotation/merged/annovar_maf_merged_filtered.qs")
# )

############################################################################
## Explore sample groups --------------
############################################################################
variant_filter_results <- qread(
    here(filter_res_dir, "variant_filter_results.qs")
)

# maf_filtered <- qread(
#   here("data/wes/annotation/merged/annovar_maf_merged_filtered.qs")
# )

maf_list <- list(
    maf = read.maf(maf = variant_filter_results$maf),
    maf_filtered = read.maf(maf = variant_filter_results$maf_filtered)
)

sample_info <- read_excel(
    here("data/clinical/DFSP-Multiomics-Sample list (updated 2024.09).xlsx")
) |> 
    select(Sample.ID, Specimen.Class, Histology.Nature)

## Sample groups
sample_groups <- list(
  u_dfsp = c("Classic", "Myxoid", "Pigmented"),
  
  pre_fst = c(
    "Pretransformed classic",
    "Pretransformed myxoid",
    "Paired classic",
    "Paired myxoid"
  ),

  post_fst = c(
    "Posttransformed FST",
    "Paired FST",
    "Paired Pleomorphic"
  ),

  fs_dfsp = c("Unpaired FST")
)

## Two samples were not appeared in the maf data
## "DFSP-139-T"    "DFSP-294-T-M1"
setdiff(
  sort(sample_info |> filter(Specimen.Class == "Tumour") |> pull(Sample.ID)),
  sort(maf_list$maf_filtered@clinical.data$Tumor_Sample_Barcode)
)

## Add the sample groups to clincial data
clinical_data <- maf_list$maf_filtered@clinical.data |> 
  left_join(
    sample_info |> 
      filter(Specimen.Class == "Tumour") |>
      rename(Tumor_Sample_Barcode = Sample.ID)
    ) |> 
  mutate(
    sample_group = case_when(
      Histology.Nature %in% sample_groups$u_dfsp ~ "u_dfsp",
      Histology.Nature %in% sample_groups$pre_fst ~ "pre_fst",
      Histology.Nature %in% sample_groups$post_fst ~ "post_fst",
      Histology.Nature %in% sample_groups$fs_dfsp ~ "fs_dfsp",
      TRUE ~ "other"
    )
  )

maf_list$maf_filtered@clinical.data <- clinical_data
maf_list$maf@clinical.data <- clinical_data

## Create the sample lists for the comparison
samples <- list(
  entire_cohort = clinical_data |> 
      pull(Tumor_Sample_Barcode),

  pre_fst_vs_post_fst = clinical_data |> 
      filter(sample_group %in% c("pre_fst", "post_fst")) |>
      pull(Tumor_Sample_Barcode),
  
  pre_fst_vs_u_dfsp = clinical_data |> 
      filter(sample_group %in% c("pre_fst", "u_dfsp")) |>
      pull(Tumor_Sample_Barcode),
  
  post_fst_vs_u_dfsp = clinical_data |>
      filter(sample_group %in% c("post_fst", "u_dfsp")) |>
      pull(Tumor_Sample_Barcode)
)

sample_group_colors <- list(
  sample_group = c(
    "u_dfsp" = "#3498db",    # Blue - for untransformed DFSP
    "pre_fst" = "#2ecc71",   # Green - for pre-transformation samples
    "post_fst" = "#e74c3c",  # Red - for post-transformation samples
    "fs_dfsp" = "#9b59b6"    # Purple - for unpaired FST
  )
)

# histology_colors <- list(
#   Histology.Nature = c(
#     # u_dfsp subtypes
#     "Classic" = "#3498db",       # Blue
#     "Myxoid" = "#2980b9",        # Darker blue
#     "Pigmented" = "#85c1e9",     # Lighter blue
    
#     # pre_fst subtypes
#     "Pretransformed classic" = "#2ecc71",  # Green
#     "Pretransformed myxoid" = "#27ae60",   # Darker green
#     "Paired classic" = "#58d68d",          # Lighter green
#     "Paired myxoid" = "#82e0aa",           # Very light green
    
#     # post_fst subtypes
#     "Posttransformed FST" = "#e74c3c",     # Red
#     "Paired FST" = "#c0392b",              # Darker red
#     "Paired Pleomorphic" = "#f1948a",      # Lighter red
    
#     # fs_dfsp subtype
#     "Unpaired FST" = "#9b59b6"             # Purple
#   )
# )

group_compare_res <- list()
top_n_genes <- 30

for (obj in names(maf_list)) {
    
    maf_obj <- maf_list[[obj]]

    for (i in names(samples)) {

        ## Get the current sample subset
        sample_list <- samples[[i]]

        maf_subset <- subsetMaf(
            maf = maf_obj,
            tsb = sample_list,
            verbose = FALSE
        )

        n_samples <- nrow(maf_subset@clinical.data)
            
        ## Perform group comparison
        if (grepl("_vs_", i)) {

            groups <- str_split(i, "_vs_")[[1]]

            group1 <- groups[1]

            group2 <- groups[2]

            ## Extract sample IDs for each group
            group1_samples <- maf_subset@clinical.data$Tumor_Sample_Barcode[
                maf_subset@clinical.data[["sample_group"]] == group1]

            group2_samples <- maf_subset@clinical.data$Tumor_Sample_Barcode[
                maf_subset@clinical.data[["sample_group"]] == group2]

            ## Check if we have enough samples in each group
            if (length(group1_samples) == 0 || length(group2_samples) == 0) {
                warning("Not enough samples in at least one group for comparison")
                return(comparison_result = NULL)
            }

            ## Perform statistical comparison
            comparison_result <- mafCompare(
                    m1 = subsetMaf(
                        maf_subset, tsb = group1_samples, verbose = FALSE
                    ),
                    m2 = subsetMaf(
                        maf_subset, tsb = group2_samples, verbose = FALSE
                    ),
                    m1Name = group1,
                    m2Name = group2,
                    minMut = 2
                )
        
        ## Extract significant results (adjust p-value threshold as needed)
        sig_res <- comparison_result$results |> 
            filter(pval < 0.01 & or > 1) |> 
            # filter(adjPval < 0.05) |>
            arrange(adjPval)
        
        message(paste("Number of significant genes:", nrow(sig_res)))

        group_compare_res[[obj]][[i]] <- comparison_result$results
        
        ## Visualize the results
        # pdf_file <- here(fig_dir, paste0("forestplot_", i ,".pdf"))
        # pdf(file = pdf_file, width = 8, height = 6)
        # forestPlot(mafCompareRes = comparison_result, pVal = 0.05)
        # dev.off()

        }
        
        # Plot the oncoplot
        MafOncoPlot(
            maf = maf_subset,
            top_n_genes = top_n_genes,
            clinical_features = "sample_group",
            annotation_colors = sample_group_colors,
            sort_by_annotation = FALSE,
            show_sample_names = FALSE,
            remove_non_mutated = FALSE,
            title = paste0(
                "Mutation Landscape: ", i, 
                " n = , ", n_samples, 
                " top ", top_n_genes, "genes)"
            ),
            font_size = 0.7,
            width = 10,
            height = 8,
            fig_dir = filter_res_dir,
            fig_name = paste0("oncoplot_", obj, "_", i),
            is_pdf = TRUE
        )
    }
    
    MafOncoPlot(
            maf = maf_object,
            top_n_genes = top_n_genes,
            clinical_features = "sample_group",
            annotation_colors = sample_group_colors,
            sort_by_annotation = FALSE,
            show_sample_names = FALSE,
            remove_non_mutated = FALSE,
            title = paste0(
                "Mutation Landscape: ", i, 
                " n = , ", n_samples, 
                " top ", top_n_genes, "genes)"
            ),
            font_size = 0.7,
            width = 10,
            height = 8,
            fig_dir = filter_res_dir,
            fig_name = paste0("oncoplot_", obj, "_", i),
            is_pdf = TRUE
    )
}

# Save the group comparison results
write_xlsx(
    group_compare_res$maf_filtered, 
    here(filter_res_dir, "post_filter_group_comparison_results.xlsx")
)

write_xlsx(
    group_compare_res$maf, 
    here(filter_res_dir, "pre_filter_group_comparison_results.xlsx")
)

