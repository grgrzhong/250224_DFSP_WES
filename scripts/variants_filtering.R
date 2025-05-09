
# Load required libraries
source(here::here("lib/R/study_lib.R"))

##############################################################################
## Merge the annovar outputs ----------------------
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

###########################################################################
## Annotate variants with cancer hotspot info --------------------
###########################################################################
## 1. Use cancer hotspot databases (e.g., COSMIC, OncoKB, or custom hotspot lists) to prioritize variants located in known cancer-associated regions.
## 2. Filter variants to retain only those that overlap with hotspot regions, as these are more likely to have clinical significance.
maf <- qread(here("data/wes/annotation/merged/annovar_maf_merged.qs")) |>
    separate(AD, into = c("RAD", "VAD"), sep = ",", remove = TRUE) |>
    mutate(
        VAD = as.numeric(VAD),
        DP = as.numeric(DP),
        AF = as.numeric(AF),
        gnomAD_exome_ALL = as.numeric(
            replace(gnomAD_exome_ALL, gnomAD_exome_ALL == ".", NA)
        )
    ) |>
    as_tibble()

# Add cancer hotspots info
maf_hotspot <- AddCancerHotspot(
    maf = maf,
    qvalue = 0.05,
    median_allele_freq_rank = 0.5,
    log10_pvalue = NULL
)

############################################################################
## Common practise for variants filtering --------------
############################################################################
# Sumarize the variants before filtering
summary_before_filter <- maf_hotspot |>
    count(Tumor_Sample_Barcode, name = "n_variants_before") |>
    summarise(
        across(
            n_variants_before,
            list(
                total = ~ sum(.x),
                mean = ~ mean(.x, na.rm = TRUE),
                median = ~ median(.x, na.rm = TRUE)
            ),
            .names = "before_filter_{.fn}"
        )
    )

message(
    sprintf("Before-filter total variants   = %d", summary_before_filter$before_filter_total), "\n",
    sprintf("Before-filter mean variants    = %.2f", summary_before_filter$before_filter_mean), "\n",
    sprintf("Before-filter median variants  = %.2f", summary_before_filter$before_filter_median), "\n"
)

## Check the columns
# SIFT_pred (D = Damaging, T = Tolerated)
# Polyphen2_HDIV_pred/Polyphen2_HVAR_pred (D = Probably damaging, P = Possibly damaging, B = Benign)
# MutationTaster_pred (A = disease causing automatic, D = disease causing, N = polymorphism, P = polymorphism automatic)
# MutationAssessor_pred (H = high, M = medium, L = low, N = neutral)
# FATHMM_pred (D = Damaging, T = Tolerated)
# PROVEAN_pred (D = Deleterious, N = Neutral)
# MetaSVM_pred/MetaLR_pred (D = Damaging, T = Tolerated)
# M-CAP_pred (D = Damaging, T = Tolerated)
# CLNSIG (ClinVar clinical significance)

colnames(maf)
## mutation type
unique(maf$Variant_Classification)
unique(maf$ExonicFunc.refGene)
unique(maf$Func.refGene)

## Functional prediction
unique(maf$SIFT_pred)
unique(maf$Polyphen2_HDIV_pred)
unique(maf$MutationTaster_pred)
unique(maf$MutationAssessor_pred)
unique(maf$FATHMM_pred)
unique(maf$PROVEAN_pred)
unique(maf$MetaSVM_pred)
unique(maf$`M-CAP_pred`)

## Use general filtering thresholds to filter out potenital sequencing errors
# or artifacts while retaining true variants
min_rd <- 8              # Minimum read depth
min_vad <- 4             # Minimum variant allele depth, 3-10
min_vaf <- 0.05          # Minimum variant allele frequency
max_pop_freq <- 0.001    # Maximum population frequency (polymorphism)

maf_filtered <- maf_hotspot |> 
    ## Filter out variants with low read depth, variant allele depth,
    dplyr::filter(
        is.na(DP) | DP >= min_rd,
        is.na(VAD) | as.numeric(VAD) >= min_vad,
        is.na(AF) | as.numeric(AF) >= min_vaf
    ) |> 
    ## Filter out variants with high population frequency
    dplyr::filter(
        is.na(gnomAD_exome_ALL) | gnomAD_exome_ALL <= max_pop_freq
    ) |> 
    arrange(gnomAD_exome_ALL) |> 
    ## Filter out silent (synonymous) variants
    dplyr::filter(
        !(Variant_Classification %in% c("Silent")),
        !(Func.refGene %in% c("synonymous_SNV")),
        !(ExonicFunc.refGene %in% c("synonymous SNV"))

    ) |> 
    ## Inlcude Filter out non-exonic variants
    dplyr::filter(
        Func.refGene %in% c(
            "exonic", "splicing", "UTR5"
        )
    )

# Sumarize the variants after filtering
summary_after_filter <- maf_filtered |>
    count(Tumor_Sample_Barcode, name = "n_variants_after") |>
    summarise(
        across(
            n_variants_after,
            list(
                total = ~ sum(.x),
                mean = ~ mean(.x, na.rm = TRUE),
                median = ~ median(.x, na.rm = TRUE)
            ),
            .names = "after_filter_{.fn}"
        )
    )

message(
    sprintf("After-filter total variants   = %d", summary_after_filter$after_filter_total), "\n",
    sprintf("After-filter mean variants    = %.2f", summary_after_filter$after_filter_mean), "\n",
    sprintf("After-filter median variants  = %.2f", summary_after_filter$after_filter_median), "\n"
)

############################################################################
## Plot the oncopot --------------
############################################################################
maf_filtered <- read.maf(maf = maf_filtered)

## Sample groups
sample_groups <- list(
    `U-DFSP` = c("Classic", "Myxoid", "Pigmented"),
    `Pre_FST` = c(
        "Pretransformed classic",
        "Pretransformed myxoid",
        "Paired classic",
        "Paired myxoid"
    ),
    `Post_FST` = c(
        "Posttransformed FST",
        "Paired FST",
        "Paired Pleomorphic"
    ),
    `FS_DFSP` = c("Unpaired FST")
)

sample_info <- LoadSampleInfo() |> 
    filter(Specimen.Class == "Tumour") |> 
    select(
        Sample.ID, Diagnosis, Specimen.Class, Specimen.Nature, Histology.Nature,
        Somatic.Status,purity,ploidy
    ) |> 
    mutate(
        sample_group = case_when(
            Histology.Nature %in% sample_groups$`U-DFSP` ~ "U-DFSP",
            Histology.Nature %in% sample_groups$`Pre_FST` ~ "Pre-FST",
            Histology.Nature %in% sample_groups$`Post_FST` ~ "Post-FST",
            Histology.Nature %in% sample_groups$`FS_DFSP` ~ "FS-DFSP",
            TRUE ~ "Other"
        )
    ) |> 
    rename(Tumor_Sample_Barcode = Sample.ID)

## Two samples were not appeared in the maf data: "DFSP-139-T", "DFSP-294-T-M1"
setdiff(
    sort(sample_info |> pull(Sample.ID)),
    sort(maf_filtered@clinical.data$Tumor_Sample_Barcode)
)

## Add the sample group to the maf data
maf_filtered@clinical.data <- maf_filtered@clinical.data |> 
    select(Tumor_Sample_Barcode) |>
    left_join(sample_info, by = "Tumor_Sample_Barcode")

annotation_colors <- list(
    sample_group = c(
        "U-DFSP"    = "#3498db", # Blue - for untransformed DFSP
        "Pre-FST"   = "#2ecc71", # Green - for pre-transformation samples
        "Post-FST"  = "#e74c3c", # Red - for post-transformation samples
        "FS-DFSP"   = "#9b59b6"  # Purple - for unpaired FST
    ),
    Specimen.Nature = c(
        "Primary"       = "royalblue", # Blue
        "Recurrence"    = "firebrick", # Green
        "Metastasis"    = "darkorange", # Red
        "Residual"      = "#9b59b6"  # Purple
    )
)

n_samples <- nrow(getSampleSummary(maf_filtered))
top_n_genes <- 30
clinical_features <- c("sample_group", "Specimen.Nature")

MafOncoPlot(
    maf = maf_filtered,
    top_n_genes = top_n_genes,
    clinicalFeatures = clinical_features,
    annotationColor = annotation_colors,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE,
    titleText = paste0(" n = ", n_samples, ", top ", top_n_genes, " genes"),
    fontSize = 0.7,
    width = 10,
    height = 8,
    fig_dir = "figures/oncoplot",
    fig_name = "oncoplot_sample_groups"
)

vc_cols <- RColorBrewer::brewer.pal(n = 8, name = "Paired")

getGeneSummary(maf_filtered)

plotmafSummary(
    maf = maf_filtered,
    color = vc_cols,
    rmOutlier = TRUE,
    addStat = "median",
    dashboard = TRUE,
    titvRaw = FALSE
)
