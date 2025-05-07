
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

############################################################################
## General variants filtering --------------
############################################################################
maf <- qread(here("data/wes/annotation/merged/annovar_maf_merged.qs")) |> 
    separate(AD, into = c("RAD", "VAD"), sep = ",", remove = TRUE) |> 
    as_tibble()

# Sumarize the variants before filtering
summary_before_filter <- maf |>
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

cat(
    sprintf("Before-filter total variants    = %d", summary_before_filter$before_filter_total), "\n",
    sprintf("Before-filter mean variants    = %.2f", summary_before_filter$before_filter_mean), "\n",
    sprintf("Before-filter median variants  = %.2f", summary_before_filter$before_filter_median), "\n"
)

colnames(maf)
unique(maf$Variant_Classification)
unique(maf$ExonicFunc.refGene)
unique(maf$Func.refGene)

maf |> 
    select(
        Variant_Classification, Func.refGene, ExonicFunc.refGene,
        VAD, DP, AF, gnomAD_exome_ALL
    )

maf <- maf |> 
    mutate(
        VAD = as.numeric(VAD),
        DP = as.numeric(DP),
        AF = as.numeric(AF),
        gnomAD_exome_ALL = as.numeric(
            replace(gnomAD_exome_ALL, gnomAD_exome_ALL == ".", NA)
        )
    )

# Use general filtering thresholds to filter out potenital sequencing errors
# or artifacts while retaining true variants
min_rd <- 8              # Minimum read depth
min_vad <- 4             # Minimum variant allele depth, 3-10
min_vaf <- 0.05          # Minimum variant allele frequency
max_pop_freq <- 0.001    # Maximum population frequency (polymorphism)

maf_filtered <- maf |> 
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

# Cancer hotspots fitlering to prioritize the most clinically relevant variants
snv_hotspots <- read_xlsx(
    here("data/clinical/hotspots_v2.xlsx"),
    sheet ="SNV-hotspots"
)

indel_hotspots <- read_xlsx(
    here("data/clinical/hotspots_v2.xlsx"),
    sheet ="INDEL-hotspots"
)
colnames(maf_filtered)
colnames(snv_hotspots)
colnames(indel_hotspots)

## Select the high confidence hotspot variants
sig_snv_hotspots <- snv_hotspots |> 
    filter(
        qvalue < 0.05,
        Median_Allele_Freq_Rank > 0.5
    ) |> 
    select(
        Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid
    ) |> 
    mutate(
        pos_aa = as.character(Amino_Acid_Position),
        ref_aa = str_extract(Reference_Amino_Acid, "^[A-Z*]"),
        var_aa = str_extract(Variant_Amino_Acid, "^[A-Z*]")
    ) |> 
    distinct() |> 
    mutate(
        snv_hotspot = paste(
            Hugo_Symbol, ref_aa, pos_aa, var_aa, sep = "_"
        )
    )

sig_indel_hotspots <- indel_hotspots |> 
    filter(
        qvalue < 0.05,
        Median_Allele_Freq_Rank > 0.5
    ) |> 
    select(
        Hugo_Symbol, Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid
    ) |> 
    mutate(
        indel_change = str_extract(Variant_Amino_Acid, "^[^:]+")
    ) |> 
    distinct() |> 
    mutate(
        indel_hotspot = paste(
            Hugo_Symbol, indel_change, sep = "_"
        )
    )

maf_filtered <- maf_filtered |> 
    # select(
    #     Tumor_Sample_Barcode, Variant_Classification, Hugo_Symbol, aaChange
    # ) |> 
    # filter(Variant_Classification == "Missense_Mutation") |> 
    mutate(
        pos_aa = str_extract(aaChange, "[0-9]+"),
        ref_aa = str_extract(aaChange, "(?<=p\\.)[A-Z*]"),
        var_aa = str_extract(aaChange, "[A-Z*]$"),
        snv_change = paste(
            Hugo_Symbol, ref_aa, pos_aa, var_aa, sep = "_"
        )
    ) |> 
    mutate(indel_change = str_remove(aaChange, "^p\\.")) |> 
    mutate(indel_change = paste(Hugo_Symbol, indel_change, sep = "_")) |> 
    mutate(
        is_snv_hotspot = if_else(
            snv_change %in% sig_snv_hotspots$snv_hotspot,
            TRUE,
            FALSE
        ),
        is_indel_hotspot = if_else(
            indel_change %in% sig_indel_hotspots$indel_hotspot,
            TRUE,
            FALSE
        )
    ) |> 
    mutate(
        is_hotspot = if_else(
        is_snv_hotspot | is_indel_hotspot,
        TRUE,
        FALSE
    )
    ) |>
    arrange(desc(is_hotspot)) |> 
    select(-pos_aa, -ref_aa, -var_aa, -snv_change, -indel_change)

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

cat(
    sprintf("After-filter total variants    = %d", summary_after_filter$after_filter_total), "\n",
    sprintf("After-filter mean variants    = %.2f", summary_after_filter$after_filter_mean), "\n",
    sprintf("After-filter median variants  = %.2f", summary_after_filter$after_filter_median), "\n"
)

############################################################################
## Plot the oncopot --------------
############################################################################
# %% Plot the oncoplot
maf_filtered <- read.maf(maf = maf_filtered)
print(maf_filtered)
