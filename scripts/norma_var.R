## Load required libraries and functions
source(here::here("bin/R/lib/study_lib.R"))

# Original format
input_files <- "/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T/DFSP-001-T.annovar.txt"
maf <- annovarToMaf(input_files)

## Clean the sample names and Clean the AD column
maf_tbl <- maf |> 
    as_tibble() |> 
    mutate(
        Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_annovar", "")
    ) |> 
    separate(AD, into = c("RAD", "VAD"), sep = ",", remove = FALSE) |>
    mutate(
        VAD = as.numeric(VAD),
        DP = as.numeric(DP),
        AF = as.numeric(AF),
        gnomAD_exome_ALL = as.numeric(
            replace(gnomAD_exome_ALL, gnomAD_exome_ALL == ".", NA)
        )
    )

# New format with both original and new columns
input_files2 <- "/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T/DFSP-001-T.annovar_test.txt"
maf2 <- annovarToMaf(input_files2)

## Clean the sample names and process new format
maf_tbl2 <- maf2 |> 
    as_tibble() |> 
    mutate(
        Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode, "_annovar", ""),
        gnomAD_exome_ALL = as.numeric(
            replace(gnomAD_exome_ALL, gnomAD_exome_ALL == ".", NA)
        )
    )

maf_tbl2 |> 
    select(
        AD, DP, AF, 
        TUMOR_AD, TUMOR_DP, TUMOR_AF, 
        NORMAL_AD, NORMAL_DP, NORMAL_AF,
    ) |>
    # Separate original AD column
    separate(AD, into = c("RAD", "VAD"), sep = ",", remove = FALSE) |>
    # Separate TUMOR_AD and NORMAL_AD columns
    separate(TUMOR_AD, into = c("TUMOR_RAD", "TUMOR_VAD"), sep = ",", remove = FALSE, fill = "right") |>
    separate(NORMAL_AD, into = c("NORMAL_RAD", "NORMAL_VAD"), sep = ",", remove = FALSE, fill = "right") |>
    mutate(
        # Convert to numeric
        VAD = as.numeric(VAD),
        DP = as.numeric(DP),
        AF = as.numeric(AF),
        TUMOR_VAD = as.numeric(TUMOR_VAD),
        TUMOR_AF = as.numeric(TUMOR_AF),
        TUMOR_DP = as.numeric(TUMOR_DP),
        NORMAL_VAD = as.numeric(NORMAL_VAD),
        NORMAL_AF = as.numeric(NORMAL_AF),
        NORMAL_DP = as.numeric(NORMAL_DP)
    )

# Compare data
print("Original format (first 5 rows):")
maf_tbl |> 
    select(AD, VAD, DP, AF, gnomAD_exome_ALL) |>
    head(5)

print("New format with both original and new columns (first 5 rows):")
maf_tbl2 |> 
    select(AD, VAD, DP, AF, TUMOR_AD, TUMOR_VAD, TUMOR_AF, TUMOR_DP, NORMAL_AD, NORMAL_VAD, NORMAL_AF, NORMAL_DP, gnomAD_exome_ALL) |>
    head(5)

# Verification: Check if original columns match new tumor columns
print("Verification - should be TRUE for all:")
verification_result <- maf_tbl2 |> 
    filter(!is.na(TUMOR_VAD), SAMPLE_TYPE == "PAIRED") |>
    summarise(
        VAD_match = all(VAD == TUMOR_VAD, na.rm = TRUE),
        AF_match = all(abs(AF - TUMOR_AF) < 0.001, na.rm = TRUE),
        DP_match = all(DP == TUMOR_DP, na.rm = TRUE),
        AD_match = all(AD == TUMOR_AD, na.rm = TRUE),
        .groups = "drop"
    )
print(verification_result)

# Filter variants with normal VAF <= 1% AND normal VAD <= 1
if ("NORMAL_AF" %in% colnames(maf_tbl2)) {
    filtered_variants <- maf_tbl2 |>
        filter(
            SAMPLE_TYPE == "PAIRED",
            !is.na(NORMAL_AF),
            !is.na(NORMAL_VAD),
            NORMAL_AF <= 0.01,
            NORMAL_VAD <= 1
        )
    
    print(paste("Total variants:", nrow(maf_tbl2)))
    print(paste("Paired variants:", sum(maf_tbl2$SAMPLE_TYPE == "PAIRED", na.rm = TRUE)))
    print(paste("Variants with normal VAF<=1% AND VAD<=1:", nrow(filtered_variants)))
    if (nrow(maf_tbl2) > 0) {
        print(paste("Percentage passing filters:", round(nrow(filtered_variants)/nrow(maf_tbl2)*100, 2), "%"))
    }
}