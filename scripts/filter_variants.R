
## Load required libraries and functions
source(here::here("bin/R/lib/study_lib.R"))

## Merge all the annovar annoated variants -------------
annovar_dir <- "data/wes/variant_calling/mutect2_with_black_repeat_filter_new"
out_dir <- "data/wes/results/merged"

# Run only once
maf_tbl <- MergeAnnovarOutput(
    annovar_dir = annovar_dir,
    is_save = TRUE,
    save_dir = out_dir
)

## Annotate variants with cancer hotspot info
maf_tbl <- LoadMergedAnnovar(out_dir)

# maf_tbl |> 
#     select(tail(names(maf_tbl), 15), gnomAD_exome_ALL) |> 
#     filter(!is_paired_normal)

## Indels/Splice sit filtering  ----------------
## VAF >= 5% and VAD >= 4 and DP >=20 in the tumor samples
filter_data1 <- maf_tbl |> 
    filter(
        is.na(tumor_AF) | tumor_AF >= 0.05, 
        is.na(tumor_DP) | tumor_DP >= 20,
        is.na(tumor_VAD) | tumor_VAD >= 4,
        is.na(gnomAD_exome_ALL) | gnomAD_exome_ALL < 0.001
    )

# filter_data1 |> select(matches("tumor|normal|AF|VAD|DP"))

## VAF <=1% and VAD<=1 in the paired normal samples
filter_data2 <- maf_tbl |> 
    filter(!is_paired_normal) |> 
    filter(
        is.na(normal_AF) | normal_AF <= 0.01, 
        is.na(normal_VAD) | normal_VAD <= 1
    )

# filter_data2 |> select(matches("tumor|normal|AF|VAD|DP"))

## Specific for tert promotor mutations
tert_gene <- filter_data1 |> 
    filter(Hugo_Symbol == "TERT") |> 
    filter(Func.refGene %in% c("exonic", "splicing", "upstream"))

## Combine the filtered data
first_inclusion <- bind_rows(
    filter_data1,
    filter_data2,
    tert_gene
) |> 
    distinct()

## All indel/splice site mutations
unique(first_inclusion$Variant_Classification)
unique(first_inclusion$ExonicFunc.refGene)
unique(first_inclusion$Variant_Type)
unique(first_inclusion$Func.refGene)

first_inclusion_splice <- first_inclusion |> 
    filter(Func.refGene %in% c("splicing"))

first_inclusion_indels <- first_inclusion |> 
    ## "SNP" "DEL" "INS" "DNP" "ONP"
    filter(Variant_Type %in% c("DEL", "INS")) |> 
    ## "Missense_Mutation" "Silent" "In_Frame_Del" "Frame_Shift_Del" NA 
    ## "Frame_Shift_Ins" "Nonsense_Mutation" "In_Frame_Ins" "Nonstop_Mutation"
    ## "Translation_Start_Site" "Unknown" "Inframe_INDEL"  
    filter(
        Variant_Classification %in% c(
            "In_Frame_Del", 
            "In_Frame_Ins", 
            "Inframe_INDEL",
            "Frame_Shift_Del", 
            "Frame_Shift_Ins", 
            "Nonsense_Mutation", 
            "Nonstop_Mutation",
            "Unknown",
            NA,
            "NA",
            "Translation_Start_Site"
        )
    )

message(
    "Number of indels/splice site mutations: ", 
    nrow(first_inclusion_indels) + nrow(first_inclusion_splice)
)

##############################################################################
## SNV filtering  ----------------
##############################################################################
## Exclude the synonymous SNVs
first_inclusion_snv <- first_inclusion |> 
    ## >=2 base substitution also include
    filter(Variant_Type  %in% c("SNP", "DNP", "ONP")) |> 
    filter(!(ExonicFunc.refGene %in% c("synonymous SNV")))

## Deleteriousness functional impact or pathogenicity predictions
## Check if variant meets at least 3 deleterious functional impact or pathogenic predictions
## First options
second_inclusion_snv1 <- first_inclusion_snv |> 
    rowwise() |>
    mutate(
        deleterious_count = sum(
            # CADD >= 20
            (!is.na(CADD_phred) && CADD_phred >= 20),
            # VEST3 criteria
            (!is.na(VEST3_score) && VEST3_score >= 0.7) || (!is.na(VEST3_rankscore) && VEST3_rankscore >= 0.9),
            # DANN >= 0.9
            (!is.na(DANN_score) && DANN_score >= 0.9),
            # SIFT prediction = "D" (Deleterious)
            (!is.na(SIFT_pred) && SIFT_pred == "D"),
            # Polyphen2 prediction = "P"/"D" (Possibly/Probably damaging)
            (!is.na(Polyphen2_HVAR_pred) && Polyphen2_HVAR_pred %in% c("P", "D")) || 
                (!is.na(Polyphen2_HDIV_pred) && Polyphen2_HDIV_pred %in% c("P", "D")),
            # MutationTaster prediction = "A"/"D" (Disease causing automatic/Disease causing)
            (!is.na(MutationTaster_pred) && MutationTaster_pred %in% c("A", "D")),
            # NA in at least one prediction algorithm qualifies
            is.na(CADD_phred) || is.na(VEST3_score) || is.na(DANN_score) || 
                is.na(SIFT_pred) || is.na(Polyphen2_HVAR_pred) || is.na(Polyphen2_HDIV_pred) || 
                is.na(MutationTaster_pred)
        )
    ) |>
    ungroup() |>
    filter(deleterious_count >= 3)

message(
    "Number of SNVs with at least 3 deleterious functional impact or pathogenic predictions: ", 
    nrow(second_inclusion_snv1)
)

## Second options
## CLNSIG, cancer hotspot, oncoKB, COSMIC
## Filter by CLNSIG, cancer hotspot, oncoKB, COSMIC
second_inclusion_snv_clnsig <- first_inclusion_snv |>
    filter(
        CLNSIG %in% c(
            "Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"
        )
    )

## COSMIC
# unique(first_inclusion_snv$cosmic70)
second_inclusion_snv_cosmic <- first_inclusion_snv |> 
    filter(!(cosmic70 %in% c(".")))

## Cancer hotspot
snv_hotspots <- LoadCancerHotspot()[["snv_hotspots"]]

dfsp_snv_hotsopts <- first_inclusion_snv |> 
    mutate(match_change = paste0(Hugo_Symbol, "_", aaChange)) |>
    mutate(
        is_hotspot = if_else(
            match_change %in% snv_hotspots$snv_hotspot,
            TRUE,
            FALSE
        )
    ) |>
    filter(is_hotspot)

second_inclusion_snv_hotspot <- dfsp_snv_hotsopts |>
    select(-match_change, -is_hotspot)

message(
    "Number of SNVs in cancer hotspots: ",
    nrow(second_inclusion_snv_hotspot)
)

## OncoKB
oncokb_gene_list <- "data/reference/OncoKB/cancerGeneList.tsv"

oncokb_gene_list <- read_tsv(here(oncokb_gene_list)) |> 
    pull(`Hugo Symbol`)

second_inclusion_snv_oncokb <- first_inclusion_snv |>
    mutate(
        is_oncokb = if_else(
            Hugo_Symbol %in% oncokb_gene_list,
            TRUE,
            FALSE
        )
    ) |> 
    filter(is_oncokb) |> 
    select(-is_oncokb)

message(
        "Number of variants with oncoKB ", 
        nrow(second_inclusion_snv_oncokb)
)

snv_oncokb_uniq <- second_inclusion_snv_oncokb |> 
    select(Hugo_Symbol, aaChange) |> 
    distinct()

message(
    "Number of unique OncoKB SNVs: ", 
    nrow(snv_oncokb_uniq)
)

# Combine all included variants
all_included_variants <- bind_rows(
    first_inclusion_splice,
    first_inclusion_indels,
    second_inclusion_snv1,
    second_inclusion_snv_clnsig,
    second_inclusion_snv_cosmic,
    second_inclusion_snv_hotspot,
    second_inclusion_snv_oncokb
) |> 
    distinct()

message(
    "Total number of filtered variants: ", 
    nrow(maf_tbl) - nrow(all_included_variants)
)

message(
    "Total number of included variants: ", 
    nrow(all_included_variants)
)

qsave(
    all_included_variants, 
    here(out_dir, "dfsp_included_variants.qs")
)

## Print out the summary of the filtered variants
message(
    "Summary of the filtered variants:\n",
    "Total variants: ", nrow(maf_tbl), "\n",
    "Indels/Splice site mutations: ", nrow(first_inclusion_indels) + nrow(first_inclusion_splice), "\n",
    "SNVs with at least 3 deleterious functional impact or pathogenic predictions: ", nrow(second_inclusion_snv1), "\n",
    "SNVs in CLNSIG: ", nrow(second_inclusion_snv_clnsig), "\n",
    "SNVs in COSMIC: ", nrow(second_inclusion_snv_cosmic), "\n",
    "SNVs in cancer hotspots: ", nrow(second_inclusion_snv_hotspot), "\n",
    "SNVs in OncoKB: ", nrow(second_inclusion_snv_oncokb), "\n",
    "Total included variants: ", nrow(all_included_variants), "\n",
    "Total filtered variants: ", nrow(maf_tbl) - nrow(all_included_variants), "\n"
)


##############################################################################
## Explore the variants ----------------
##############################################################################
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