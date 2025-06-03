library(tidyverse)
library(glue)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description = "Create Panel of Normals for FACETS CNV analysis")
parser$add_argument("--normal_bam_dir", required = TRUE, 
                   help = "Directory containing normal BAM files")
parser$add_argument("--output_dir", required = TRUE,
                   help = "Output directory for PON files")
parser$add_argument("--snp_vcf", required = TRUE,
                   help = "dbSNP VCF file for SNP positions")
parser$add_argument("--mutect2_pon", default = NULL,
                   help = "Optional: Mutect2 PON VCF file")
parser$add_argument("--max_normals", type = "integer", default = 30,
                   help = "Maximum number of normal samples to use")
parser$add_argument("--min_normals_for_pon", type = "integer", default = 10,
                   help = "Minimum normals needed to create PON")

args <- parser$parse_args()

# Configuration
normal_bam_dir <- args$normal_bam_dir
output_dir <- args$output_dir
snp_vcf <- args$snp_vcf
mutect2_pon <- args$mutect2_pon
max_normals <- args$max_normals
min_normals_for_pon <- args$min_normals_for_pon

# Create output directory
fs::dir_create(output_dir, recurse = TRUE)

# Get normal sample lists
normal_samples <- list.files(normal_bam_dir, pattern = "\\.bam$", full.names = FALSE) %>%
    str_remove("\\.bam$")

normal_bam_paths <- fs::path(normal_bam_dir, glue("{normal_samples}.bam"))

cat(glue("Found {length(normal_samples)} normal samples\n"))

if (length(normal_samples) < min_normals_for_pon) {
    stop(glue("Need at least {min_normals_for_pon} normal samples, found {length(normal_samples)}"))
}

# Function to create FACETS PON from normal BAMs
create_facets_pon <- function(normal_bams, pon_output_file, snp_file, use_mutect2 = FALSE) {
    cat(glue("Creating FACETS PON from {length(normal_bams)} normal samples\n"))
    
    # Select subset of normals for PON
    selected_normals <- if (length(normal_bams) > max_normals) {
        sample(normal_bams, max_normals)
    } else {
        normal_bams
    }
    
    cat(glue("Using {length(selected_normals)} normal samples for PON\n"))
    
    if (use_mutect2 && !is.null(mutect2_pon)) {
        # Create temporary SNP file from Mutect2 PON
        temp_snp_file <- tempfile(fileext = ".vcf")
        
        cat("##fileformat=VCFv4.2\n", file = temp_snp_file)
        cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = temp_snp_file, append = TRUE)
        
        cmd_extract <- glue("zcat {mutect2_pon} | grep -v '^#' | awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,\".\",$4,$5,\".\",\"PASS\",\".\"}}' >> {temp_snp_file}")
        system(cmd_extract)
        
        snp_file <- temp_snp_file
        cat("Using Mutect2 PON sites for SNP positions\n")
    }
    
    # Create pooled normal using first two samples as "fake" tumor-normal pair
    if (length(selected_normals) >= 2) {
        fake_tumor <- selected_normals[1]
        fake_normals <- selected_normals[-1]
        
        normal_bam_string <- str_c(fake_normals, collapse = " ")
        
        cmd <- glue("snp-pileup -g -q15 -Q20 -P100 -r25,0 {snp_file} {pon_output_file} {normal_bam_string} {fake_tumor}")
        
        cat(glue("Running: {cmd}\n"))
        result <- system(cmd)
        
        # Clean up temporary file
        if (exists("temp_snp_file") && fs::file_exists(temp_snp_file)) {
            fs::file_delete(temp_snp_file)
        }
        
        return(result == 0 && fs::file_exists(pon_output_file))
    }
    
    return(FALSE)
}

# Create different PON variants
pon_variants <- list(
    dbsnp = list(
        file = fs::path(output_dir, "facets_pon_dbsnp.csv.gz"),
        snp_file = snp_vcf,
        use_mutect2 = FALSE
    )
)

# Add Mutect2 PON variant if available
if (!is.null(mutect2_pon) && fs::file_exists(mutect2_pon)) {
    pon_variants$mutect2 <- list(
        file = fs::path(output_dir, "facets_pon_mutect2.csv.gz"),
        snp_file = mutect2_pon,
        use_mutect2 = TRUE
    )
}

# Create PON files
pon_results <- map_dfr(names(pon_variants), ~ {
    variant_name <- .x
    variant_info <- pon_variants[[.x]]
    
    cat(glue("\n=== Creating {variant_name} PON ===\n"))
    
    success <- create_facets_pon(
        normal_bams = normal_bam_paths,
        pon_output_file = variant_info$file,
        snp_file = variant_info$snp_file,
        use_mutect2 = variant_info$use_mutect2
    )
    
    tibble(
        pon_type = variant_name,
        pon_file = variant_info$file,
        success = success,
        num_normals = length(normal_bam_paths),
        created_date = Sys.time()
    )
})

# Save PON manifest
pon_results %>%
    write_tsv(fs::path(output_dir, "pon_manifest.txt"))

# Save sample list used for PON
tibble(
    sample_id = str_remove(basename(normal_bam_paths), "\\.bam$"),
    bam_path = normal_bam_paths
) %>%
write_tsv(fs::path(output_dir, "pon_samples.txt"))

cat(glue("\nPON creation completed. Results saved in: {output_dir}\n"))
cat("Created PON files:\n")
pon_results %>%
    filter(success) %>%
    pull(pon_file) %>%
    walk(~ cat(glue("  - {.x}\n")))
