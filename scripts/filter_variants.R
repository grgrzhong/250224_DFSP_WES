
## Load required libraries and functions
source(here::here("bin/R/lib/study_lib.R"))

## Merge all the annovar annoated variants
annovar_dir <- "data/wes/variant_calling/mutect2_with_black_repeat_filter_new"

# Run only once
maf_tbl <- MergeAnnovarOutput(
    annovar_dir = annovar_dir,
    is_save = TRUE,
    save_dir = "data/wes/variant_calling/merged"
)

## Annotate variants with cancer hotspot info
maf_tbl <- LoadMergedAnnovar()

## 1. Use cancer hotspot databases (e.g., COSMIC, OncoKB, or custom hotspot lists) to prioritize variants located in known cancer-associated regions.
## 2. Filter variants to retain only those that overlap with hotspot regions, as these are more likely to have clinical significance.
maf_tbl <- AddCancerHotspot(
    maf = maf_tbl,
    qvalue = 0.05,
    median_allele_freq_rank = 0.5,
    log10_pvalue = NULL
)

## Filter the variants based on common criteria
filter_params <- list(
    # Minimum read depth
    DP = list(op = ">=", value = 8),

    # Minimum variant allele depth
    VAD = list(op = ">=", value = 4),

    # Minimum variant allele frequency
    AF = list(op = ">=", value = 0.05),

    # Maximum population frequency
    gnomAD_exome_ALL = list(op = "<=", value = 0.01),

    # Variant classifications mutation type
    Variant_Classification = list(
        op = "out", 
        value = c("Silent")
        # other options are included: 
        # NA, "Unknown", "Missense_Mutation", "In_Frame_Del,"Inframe_INDEL",
        # "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation",
        # "In_Frame_Ins", "Translation_Start_Site", "Nonstop_Mutation",
    ),

    # Variant classifications exonic function
    ExonicFunc.refGene = list(
        op = "out", 
        value = c("synonymous_SNV")
        # other options are included: 
        # ".", "nonsynonymous SNV", "nonframeshift deletion",
        # "frameshift deletion", "frameshift insertion", "stopgain", 
        # "nonframeshift insertion","startloss", "stoploss", "unknown", 
        # "nonframeshift substitution"
    ),

    # Variant classifications pathogenicity predictions
    # Which tools to be used?

    # Variant classifications functional consequences
    Func.refGene = list(
        op = "in", 
        value = c("exonic", "UTR5", "splicing", "upstream")
        # other options are excluded:
        # "intronic", "ncRNA_exonic", "splicing", "ncRNA_splicing",
        # "ncRNA_intronic", UTR3", "upstream", "intergenic", "downstream", 
    )
)


# availabel prediction tools for functional consequences
maf_tbl |> select(matches("pred")) |> colnames()


maf_filter <- FilterMergedMaf(
    maf_tbl = maf_tbl,
    filter_params = filter_params
)

# maf_filter |> 
    # filter(Variant_Classification == "Silent")
    # filter(Func.refGene == "ncRNA_exonic")
    # filter(ExonicFunc.refGene == "synonymous_SNV")

## Find the enriched variants in FST
maf_obj <- read.maf(maf = maf_filter)

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
# setdiff(
#     sort(sample_info |> pull(Tumor_Sample_Barcode)),
#     sort(maf_obj@clinical.data$Tumor_Sample_Barcode)
# )

## Add the sample group to the maf data
maf_obj@clinical.data <- maf_obj@clinical.data |> 
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
        "Primary"       = "#ff9800", # Orange
        "Recurrence"    = "#009688", # Teal
        "Metastasis"    = "#795548", # Brown
        "Residual"      = "#607d8b"  # Blue-gray
    )
)

n_samples <- nrow(getSampleSummary(maf_obj))
top_n_genes <- 30
clinical_features <- c("sample_group", "Specimen.Nature")

oncoplot(
    maf = maf_obj,
    clinicalFeatures = clinical_features,
    annotationColor = annotation_colors,
    sortByAnnotation = TRUE,
    showTumorSampleBarcodes = FALSE,
    removeNonMutated = FALSE
)

# MafOncoPlot(
#     maf = maf_obj,
#     top_n_genes = top_n_genes,
#     clinicalFeatures = clinical_features,
#     annotationColor = annotation_colors,
#     sortByAnnotation = TRUE,
#     showTumorSampleBarcodes = FALSE,
#     removeNonMutated = FALSE,
#     titleText = paste0(" n = ", n_samples, ", top ", top_n_genes, " genes"),
#     fontSize = 0.7,
#     width = 10,
#     height = 8,
#     fig_dir = "figures/oncoplot",
#     fig_name = "oncoplot_sample_groups"
# )