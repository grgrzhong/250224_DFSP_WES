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

# tumour samples without matched normal
tumour_samples <- vcf_samplesheet |> 
    filter(status == 1) |> 
    select(sample_id, vcf, tbi, f1r2)

message(length(assume_normal_ids))

normal_ids <- grep(
    "N",
    vcf_samplesheet$sample_id, 
)

sample_info <- read_csv(here("data/wes/sample_info/samplesheet.csv")) |> 
    select(patient, sample, status)

assume_normal_ids <- paste0(sample_info$patient, "-N") |> unique()
avail_normal_ids <- sample_info |> 
    filter(status == 0) |> 
    pull(sample)

no_normal_patient <- setdiff(
    assume_normal_ids,
    avail_normal_ids
)
no_normal_patient <- str_replace_all(
    no_normal_patient,
    "-N",
    ""
)

tumour_only_samples <- sample_info |> 
    filter(patient %in% no_normal_patient) |> 
    pull(sample)
    
# Write tumour-only samples to a text file
write_lines(
    tumour_only_samples,
    here("data/wes/sample_info/tumour_only_samples.txt")
)
