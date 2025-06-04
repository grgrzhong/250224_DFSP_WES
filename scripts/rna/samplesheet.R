library(tidyverse)
library(here)

## Output directory for the sample sheet
out_dir <- here("data/rna/sample_info")

## Create the sample sheet for raw fastq files ---------------------------
fastq_dir <- here("/mnt/m/RNA-seq/STUMP/primary_seq")
fastq_paths <- dir_ls(fastq_dir, glob = "*.fastq.gz", recurse = TRUE)
fastq_files <- fastq_paths |> 
    path_file() |> 
    path_ext_remove()

# all_samples <- sapply(
#     fastq_files,
#     function(x) str_split(x, "_")[[1]][1]
# ) |> unique()

fastq_samplesheet <- tibble(
    files = fastq_files,
    paths = fastq_paths
) |> 
    separate(
        files, 
        into = c("sample", "read"), 
        sep = "_", 
        extra = "merge", 
        fill = "right"
    ) |>
    mutate(
        read = case_when(
            str_detect(read, "1") ~ "1", 
            str_detect(read, "2") ~ "2", 
            TRUE ~ "unknown"
        )
    ) |> 
    pivot_wider(
        names_from = read, 
        values_from = paths, 
        names_prefix = "fastq_"
    )

write_csv(
    fastq_samplesheet, 
    file = here(out_dir, "fastq_samplesheet.csv")
)

## Create the sample sheet for trimmed fastq files ----------------
fastq_trimmed_dir <- here("/mnt/m/RNA-seq/STUMP/Input-trimmed")

fastq_trimmed_paths <- dir_ls(fastq_trimmed_dir, glob = "*.fastq.gz", recurse = TRUE)
fastq_trimmed_files <- fastq_trimmed_paths |> 
    path_file() |> 
    path_ext_remove()

fastq_trimmed_samplesheet <- tibble(
    files = fastq_trimmed_files,
    paths = fastq_trimmed_paths
) |> 
    separate(
        files, 
        into = c("sample", "read"), 
        sep = "_", 
        extra = "merge", 
        fill = "right"
    ) |>
    mutate(
        read = case_when(
            str_detect(read, "R1") ~ "1", 
            str_detect(read, "R2") ~ "2", 
            TRUE ~ "unknown"
        )
    ) |> 
    pivot_wider(
        names_from = read, 
        values_from = paths, 
        names_prefix = "fastq_trimmed_"
    )

write_csv(
    fastq_trimmed_samplesheet, 
    file = here(out_dir, "fastq_trimmed_samplesheet.csv")
)

## Save all the samples in a single samplesheet ----------------
samplesheet <- fastq_samplesheet |> 
    left_join(
        fastq_trimmed_samplesheet, 
        by = "sample"
    )

write_csv(
    samplesheet, 
    file = here(out_dir, "samplesheet.csv")
)

write_csv(
    samplesheet |> filter(sample=="S56"), 
    file = here(out_dir, "test.csv")
)

## Check for samples not found in trimmed fastq files
not_trmmed_samples <- setdiff(fastq_samplesheet$sample, fastq_trimmed_samplesheet$sample) |> 
    str_c(collapse = ", ")

message("Samples not found in trimmed fastq files: ", not_trmmed_samples)
