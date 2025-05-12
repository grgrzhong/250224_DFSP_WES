# Load required libraries
source(here::here("lib/R/study_lib.R"))

vcf_dir <- here("data/wes/variant_calling/mutect2")
bam_dir <- here("data/wes/preprocessing/recalibrated")

# Find the recalibrated BAM files
bam_files <- dir_ls(bam_dir, recurse = TRUE, glob = "*recalibrated.bam")
bai_files <- dir_ls(bam_dir, recurse = TRUE, glob = "*recalibrated.bai")

# Find the unfiltered VCF files
vcf_files <- dir_ls(vcf_dir, recurse = TRUE, glob = "*unfiltered.vcf.gz")
idx_files <- dir_ls(vcf_dir, recurse = TRUE, glob = "*unfiltered.vcf.gz.tbi")
f1r2_files <- dir_ls(vcf_dir, recurse = TRUE, glob = "*f1r2.tar.gz")

bam_sample_ids <- bam_files |> 
    str_extract("(?<=/)[^/]+(?=_recalibrated.bam)") |>
    str_replace_all("_", "-")

bam_patient_ids <- bam_sample_ids |> 
    str_extract("^[^-]+-[^-]+")

bam_patient_status <- if_else(
    grepl("N", bam_sample_ids, ignore.case = TRUE),
    0, 1 
)

vcf_samples_ids <- vcf_files |>
    str_extract("(?<=/)[^/]+(?=_unfiltered.vcf.gz)") |>
    str_replace_all("_", "-")

bam_tbl <- tibble(
    patient_id = bam_patient_ids,
    sample_id = bam_sample_ids,
    status = bam_patient_status,
    bam = as.character(bam_files),
    bai = as.character(bai_files)
)

vcf_tbl <- tibble(
    sample_id = vcf_samples_ids,
    vcf = as.character(vcf_files),
    tbi = as.character(idx_files),
    f1r2 = as.character(f1r2_files)
)

samplesheet <- bam_tbl |> 
    left_join(vcf_tbl)

vcf_samplesheet <- samplesheet |> 
    filter(!is.na(vcf))

write_excel_csv(
    vcf_samplesheet,
    here("data/wes/csv/mutect2_call_vcf.csv")
)
