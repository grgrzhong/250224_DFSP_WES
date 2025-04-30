library(maftools)
library(qs)
library(fs)
library(here)
library(tidyverse)
library(readxl)

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
maf <- qread(here("data/wes/annotation/merged/annovar_maf_merged.qs"))

maf <- maf |> 
  separate(AD, into = c("RAD", "VAD"), sep = ",", remove = TRUE) |> 
  as_tibble()

colnames(maf)
unique(maf$Variant_Classification)
unique(maf$Func.refGene)
unique(maf$ExonicFunc.refGene)

# Define filtering thresholds to filter out potenital sequencing errors or 
# artifacts while retaining true variants
filter_params <- list(
  min_depth = 10,                   # Minimum read depth
  min_vaf = 0.05,                   # Minimum variant allele frequency
  min_vad = 5,                      # Minimum variant allele depth, 3-10
  max_population_freq = 0.001,       # Maximum population frequency
  exclude_classifications = c(      # Variant classifications to exclude
    "Silent", "Intron", "3'UTR", 
    "5'UTR", "3'Flank", "5'Flank", "IGR"
  ),
  exclude_exonic_func = c(          # Exonic functions to exclude
    "synonymous_SNV"
  ),
  include_func = c(                 # functions to include
    "exonic", "splicing"
  )
)

# filter_params1, min_depth = 20 is too stringent
filter_params1 <- filter_params 
filter_params1$min_depth <- 20

filter_res <- TestVariantFilter(
  filter_params = filter_params1,
  maf,
  is_save = TRUE,
  fig_dir = "figures/variant_filtering",
  fig_name = "filter_params1"
)

# filter_params2
filter_params2 <- filter_params 
filter_params2$min_depth <- 10

filter_res <- TestVariantFilter(
  filter_params = filter_params2,
  maf,
  is_save = TRUE,
  fig_dir = "figures/variant_filtering",
  fig_name = "filter_params2"
)

# filter_params3
filter_params3 <- filter_params 
filter_params3$min_depth <- 10
filter_params3$min_vad <- 3

filter_res <- TestVariantFilter(
  filter_params = filter_params3,
  maf,
  is_save = TRUE,
  fig_dir = "figures/variant_filtering",
  fig_name = "filter_params3"
)

# filter_params4
filter_params4 <- filter_params 
filter_params4$min_depth <- 10
filter_params4$min_vad <- 3
filter_params4$min_vaf <- 0.04

filter_res <- TestVariantFilter(
  filter_params = filter_params4,
  maf,
  is_save = TRUE,
  fig_dir = "figures/variant_filtering",
  fig_name = "filter_params4"
)

# Create a MAF object from filtered data for further analysis
maf_filtered <- read.maf(maf = filter_res$maf_filtered)

# Save filtered data
qsave(
  maf_filtered, 
  here("data/wes/annotation/merged/annovar_maf_merged_filtered.qs")
)

############################################################################
## Explore sample groups --------------
############################################################################
maf_filtered <- qread(
  here("data/wes/annotation/merged/annovar_maf_merged_filtered.qs")
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
  sort(maf_filtered@clinical.data$Tumor_Sample_Barcode)
)

## Add the sample groups to clincial data
clinical_data <- maf_filtered@clinical.data |> 
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

maf_filtered@clinical.data <- clinical_data

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

histology_colors <- list(
  Histology.Nature = c(
    # u_dfsp subtypes
    "Classic" = "#3498db",       # Blue
    "Myxoid" = "#2980b9",        # Darker blue
    "Pigmented" = "#85c1e9",     # Lighter blue
    
    # pre_fst subtypes
    "Pretransformed classic" = "#2ecc71",  # Green
    "Pretransformed myxoid" = "#27ae60",   # Darker green
    "Paired classic" = "#58d68d",          # Lighter green
    "Paired myxoid" = "#82e0aa",           # Very light green
    
    # post_fst subtypes
    "Posttransformed FST" = "#e74c3c",     # Red
    "Paired FST" = "#c0392b",              # Darker red
    "Paired Pleomorphic" = "#f1948a",      # Lighter red
    
    # fs_dfsp subtype
    "Unpaired FST" = "#9b59b6"             # Purple
  )
)

fig_dir <- "figures/oncoplot"

for (i in names(samples)) {

    ## Get the current sample subset
    sample_list <- samples[[i]]

    maf_subset <- subsetMaf(
      maf = maf_filtered,
      tsb = sample_list,
      verbose = FALSE
    )

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
              m1 = subsetMaf(maf_subset, tsb = group1_samples, verbose = FALSE),
              m2 = subsetMaf(maf_subset, tsb = group2_samples, verbose = FALSE),
              m1Name = group1,
              m2Name = group2,
              minMut = 2
          )
      
      ## Extract significant results (adjust p-value threshold as needed)
      sig_res <- comparison_result$results |> 
          filter(pval < 0.01 & or > 1) |> 
          arrange(adjPval)
      
      message(paste("Number of significant genes:", nrow(sig_res)))

      ## Save the results
      res_dir <- here("results/maf_group_comparsion")
      fs::dir_create(res_dir)
      write_excel_csv(
          comparison_result$results, 
          file = here(res_dir, paste0("group_comparison_", i,".csv"))
      )
      
      message(paste("Saved comparison results to:", here(res_dir)))

      ## Visualize the results
      pdf_file <- here(fig_dir, paste0("forestplot_", i ,".pdf"))
      pdf(file = pdf_file, width = 8, height = 6)
      forestPlot(mafCompareRes = comparison_result, pVal = 0.05)
      dev.off()

    }
    
    # Plot the oncoplot
    p <- MafOncoPlot(
      maf = maf_subset,
      top_n_genes = 30,
      clinical_features = "sample_group",
      annotation_colors = sample_group_colors,
      sort_by_annotation = FALSE,
      show_sample_names = FALSE,
      remove_non_mutated = FALSE,
      title = paste0("Mutation Landscape: ", i),
      font_size = 0.7,
      width = 10,
      height = 8,
      fig_dir = fig_dir,
      fig_name = paste0("oncoplot_", i)
    )

}
