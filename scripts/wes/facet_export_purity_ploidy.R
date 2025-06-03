
## Load required libraries
source(here::here("bin/R/lib/study_lib.R"))

## List all cnv facet VCF files
cnv_facet_dir <- here("data/wes/variant_calling/cnv/facets")

vcf_files <- list.files(
    path = cnv_facet_dir,
    pattern = "\\.vcf\\.gz$",
    recursive = TRUE,
    full.names = TRUE
)
message(paste("Found", length(vcf_files), "VCF files"))

## Function to extract values from VCF header
extract_facet_info <- function(vcf_file) {
    
    header <- system(
        paste("zcat", shQuote(vcf_file), "| grep '^##'"), intern = TRUE
    )
    
    get_val <- function(key) {
        val <- sub(
            paste0("##", key, "="), "", header[grep(paste0("^##", key, "="), header)]
        )
        if (length(val) == 0) {
            return(NA)
        }
        val
    }
    
    sample <- basename(dirname(vcf_file))
    
    tibble(
        sample = sample,
        purity = as.numeric(get_val("purity")),
        ploidy = as.numeric(get_val("ploidy")),
        dipLogR = as.numeric(get_val("dipLogR")),
        est_insert_size = as.numeric(get_val("est_insert_size")),
        emflags = get_val("emflags")
    )
}

## Extract info for all files
facet_info <- list_rbind(lapply(vcf_files, extract_facet_info))

write_excel_csv(
    facet_info, 
    file = here(cnv_facet_dir, "cnv_facet_purity_ploidy.csv"),
)

## Add the ploidy and purity values to the sample sheet
sample_info <- read_excel(
    here("data/clinical/DFSP-Multiomics-Sample list (updated 2024.09).xlsx")
)

sample_info <- sample_info |> 
    left_join(
        facet_info |> select(sample, purity, ploidy),
        by = c("Sample.ID" = "sample")
    )

write_xlsx(
    sample_info,
    here("data/clinical/DFSP_multiomics_sample_list_updated_20250509.xlsx")
)
