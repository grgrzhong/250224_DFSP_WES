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

maf <- annovarToMaf(input_files)

## Clean the sample names
maf <- maf |> 
  mutate(
    Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_annovar", "")
  )

qsave(maf, here("data/wes/annotation/merged/merged_annovar_maf.qs"))

##############################################################################
## Merge the funcotator output --------------
##############################################################################
# input_dir <- here("data/wes/annotation/funcotator")
# input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*annotated.tsv")

# input_data <- map2(
#   input_files, names(input_files),
#   ~ {
#     maf_data <- read_tsv(.x, show_col_types = FALSE)
#     maf_data$Tumor_Sample_Barcode <- basename(.y) |> str_remove("\\..*$")

#     return(maf_data)
#   }
# )

# maf <- merge_mafs(mafs = input_data, verbose = TRUE)

# qsave(maf, here("data/wes/annotation/merged/merged_funcotator_maf.qs"))

##############################################################################
## Filtering the variants --------------
##############################################################################
maf <- qread(here("data/wes/annotation/merged/merged_annovar_maf.qs"))

maf <- maf |> 
  separate(AD, into = c("RAD", "VAD"), sep = ",", remove = TRUE) |> 
  as_tibble()

# Define filtering thresholds
filter_params <- list(
  min_depth = 20,                   # Minimum read depth
  min_vaf = 0.05,                   # Minimum variant allele frequency
  max_population_freq = 0.01,       # Maximum population frequency
  exclude_classifications = c(      # Variant classifications to exclude
    "Silent", "Intron", "3'UTR", 
    "5'UTR", "3'Flank", "5'Flank", "IGR"
  ),
  exclude_exonic_func = c(          # Exonic functions to exclude
    "synonymous SNV", "unknown"
  )
)

# Test the filtering
test_fiter <- TestVariantFilter(
  filter_params,
  maf,
  is_save = TRUE,
  fig_dir = "figures/annovar",
  fig_name = "compare_filter"
)

# Create a MAF object from filtered data for further analysis
maf_filtered_obj <- read.maf(maf = test_fiter$maf_filtered)

# Save filtered data
qsave(maf_filtered, here("data/wes/annotation/filtered/filtered_maf.qs"))

##############################################################################
## Explore sample groups --------------
##############################################################################
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

## Do variants filtering
maf_filtered <- maf |> 
  as_tibble()

maf |> colnames()

# getSampleSummary(maf_data)
# getGeneSummary(maf_data)
# getClinicalData(maf_data)
# getFields(maf_data)
# write.mafSummary(maf = maf_data, basename = "maf_data")

plotmafSummary(
  maf = maf_data,
  rmOutlier = TRUE,
  addStat = "median",
  dashboard = TRUE,
  titvRaw = FALSE
)

oncoplot(maf = maf_data, top = 20)

maf_data.titv <- titv(maf = maf_data, plot = FALSE, useSyn = TRUE)
# plot titv summary
plotTiTv(res = maf_data.titv)

maf_data@data |>
  as_tibble() |>
  arrange(Hugo_Symbol) |>
  filter(
    Hugo_Symbol %in% c("MYOD1", "SUZ12", "NF1")
  )

lollipopPlot(
  maf = maf_data,
  gene = "MYOD1",
  AACol = "Protein_Change",
  showMutationRate = TRUE,
  labelPos = "all"
)

lollipopPlot(
  maf = maf_data,
  gene = "SUZ12",
  AACol = "Protein_Change",
  showMutationRate = TRUE,
  labelPos = "all"
)

maf_data@data |>
  as_tibble() |>
  arrange(Hugo_Symbol) |>
  filter(
    # Variant_Classification %in% c("Splice_Site"),
    Hugo_Symbol %in% c("MYOD1", "SUZ12", "NF1")
  ) |>
  # view()
  select(
    Hugo_Symbol, Tumor_Sample_Barcode,
    Genome_Change, cDNA_Change, Protein_Change, Codon_Change
  )
