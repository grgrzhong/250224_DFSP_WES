## Load required libraries
source(here::here("lib/R/study_lib.R"))

#############################################################################
## Export bam csv ----
#############################################################################
## Standard csv columns
bam_std_cols <- c(
    "patient", "sample", "status", "bam", "bai"
)

## Get the include samples from the bam files
# bam_dir <- here("/mnt/m/WES/DFSP/BAM")
bam_dir <- here("data/wes/preprocessing/recalibrated")
bam_files <- dir_ls(bam_dir, recurse = TRUE, glob = "*recalibrated.bam")
samples <- path_file(bam_files) |>
    path_ext_remove() |> 
    str_remove("_recalibrated")

message("Found ", length(samples), " samples")
bam_tbl <- tibble(
    sample = samples,
    patient = str_extract(samples, "DFSP-\\d+"),
    status = if_else(str_detect(samples, "-T"), 1, 0),
    bam = bam_files,
    bai = str_replace(bam_files, ".bam", ".bai")
) |> 
    select(all_of(bam_std_cols))

## write to CSV
out_dir <- "data/wes/sample_info"

write_excel_csv(bam_tbl, here(out_dir, "bam_all_samples.csv"))

#############################################################################
## Export fastq csv ----
#############################################################################
fastq_std_cols <- c(
    "patient", "sample", "status", "fastq_1", "fastq_2"
)

## Get the matched fastq files
fastq_dir <- here("/mnt/m/WES/DFSP/Raw")
fastq_files <- dir_ls(fastq_dir, recurse = TRUE, glob = "*fastq.gz")
fastq_names <- fastq_files |> 
    path_file() |> 
    path_ext_remove() |> 
    path_ext_remove()  

# Create dataframe with file info
fastq_tbl <- tibble(
    filepath = fastq_files,
    filename = fastq_names
) |>
    # Split into R1/R2
    mutate(read = str_extract(filename, "_[12]")) |> 
    mutate(sample = str_remove(filename, "_[12]")) |> 
    mutate(read = str_remove(read, "_"))

message("Found ", nrow(fastq_tbl), " fastq files")
message("Found ", length(unique(fastq_tbl$sample)), " samples")

## Check the fastq files
check_data <- fastq_tbl |> 
    group_by(sample) |>
    summarise(
        n_fastq = n(),
        n_fastq_1 = sum(read == "1"),
        n_fastq_2 = sum(read == "2")
    ) |>
    mutate(
        status = case_when(
            n_fastq < 2 ~ "missing files",
            n_fastq > 2 ~ "extra files",
            TRUE ~ "ok"
        )
    ) |>
    filter(status != "ok")

fastq_tbl <- fastq_tbl |> 
    pivot_wider(
        id_cols = sample,
        names_from = read,
        values_from = filepath,
        names_prefix = "fastq_"
    ) |>
    # Add patient ID and status
    mutate(
        patient = str_extract(sample, "DFSP-\\d+"),
        status = if_else(str_detect(sample, "-T"), 1, 0)
    ) |>
    # Reorder columns
    select(all_of(fastq_std_cols))

## Check if any not matched with bam files
table(samples %in% fastq_csv$sample)

# Write to csv
write_csv(
    fastq_tbl, 
    here(out_dir, "fastq_all_samples.csv")
)

write_csv(
    fastq_tbl |> filter(sample == "DFSP-001-T"),
    here(out_dir, "fastq_test1.csv")
)
