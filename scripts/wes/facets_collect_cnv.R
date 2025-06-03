## Load required libraries
source(here::here("bin/R/lib/study_lib.R"))

## List all CNV facet annotated TSV files
facets_dir1 <- here("data/wes/variant_calling/cnv/facets")
facets_dir2 <- here("data/wes/variant_calling/cnv/facets_test_normal")

facets_files1 <- dir_ls(facets_dir1, recurse = TRUE, glob = "*.annotated.tsv")

facets_files2 <- dir_ls(facets_dir2, recurse = TRUE, glob = "*.annotated.tsv")

facets_files <- c(facets_files1, facets_files2)

message("Found ", length(facets_files), " CNV facet annotated TSV files.")

## Load all CNV facet annotated TSV files
facets_data <- map(
    facets_files, 
    ~ read_tsv(
        .x, 
        col_types = cols(.default = "c"), 
        show_col_types = FALSE
    ),
    .progress = TRUE
)

facets_tbl <- list_rbind(facets_data)

message("Found ", length(unique(facets_tbl$Sample.ID)), " samples")

## Save the combined data
write_csv(
    facets_tbl,
    here(facets_dir, "cnv_facets_wes_all_tumours.csv")
)

## Prepare the input the maftools
cnv_facets <- facets_data |> 
    select(
        Sample.ID, 
        Gene,
        SVTYPE,
    )