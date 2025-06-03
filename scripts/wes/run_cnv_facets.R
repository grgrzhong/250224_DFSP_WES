library(facets)
library(pctGCdata)
library(tidyverse)
library(parallel)
library(glue)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description = "Run FACETS CNV analysis using matched normals or PON")
parser$add_argument("--tumor_bam_dir", required = TRUE,
                   help = "Directory containing tumor BAM files")
parser$add_argument("--output_dir", required = TRUE,
                   help = "Output directory for FACETS results")
parser$add_argument("--snp_vcf", required = TRUE,
                   help = "dbSNP VCF file for SNP positions")

# Mutually exclusive group for normal source
normal_group <- parser$add_mutually_exclusive_group(required = TRUE)
normal_group$add_argument("--pon_file", 
                         help = "Pre-generated PON file from facets_create_pon.R")
normal_group$add_argument("--normal_bam_dir",
                         help = "Directory containing matched normal BAM files")

parser$add_argument("--sample_sheet", 
                   help = "CSV file with tumor_id,normal_id columns (required if using --normal_bam_dir)")
parser$add_argument("--genome_build", default = "hg38",
                   help = "Genome build (hg19 or hg38)")
parser$add_argument("--cval", type = "integer", default = 150,
                   help = "Critical value for segmentation")
parser$add_argument("--parallel", type = "integer", default = 1,
                   help = "Number of parallel processes")

args <- parser$parse_args()

# Configuration
tumor_bam_dir <- args$tumor_bam_dir
output_dir <- args$output_dir
snp_vcf <- args$snp_vcf
pon_file <- args$pon_file
normal_bam_dir <- args$normal_bam_dir
sample_sheet <- args$sample_sheet
genome_build <- args$genome_build
cval <- args$cval
n_cores <- args$parallel

# Validate inputs
if (!fs::dir_exists(tumor_bam_dir)) {
    stop(glue("Tumor BAM directory not found: {tumor_bam_dir}"))
}

# Validate normal source
use_pon <- !is.null(pon_file)
use_matched_normals <- !is.null(normal_bam_dir)

if (use_pon) {
    if (!fs::file_exists(pon_file)) {
        stop(glue("PON file not found: {pon_file}"))
    }
    cat(glue("Using PON approach with: {pon_file}\n"))
} else if (use_matched_normals) {
    if (!fs::dir_exists(normal_bam_dir)) {
        stop(glue("Normal BAM directory not found: {normal_bam_dir}"))
    }
    if (is.null(sample_sheet) || !fs::file_exists(sample_sheet)) {
        stop("Sample sheet required when using matched normals")
    }
    cat(glue("Using matched normals from: {normal_bam_dir}\n"))
}

# Create output directory
fs::dir_create(output_dir, recurse = TRUE)

# Get tumor sample lists
tumor_samples <- list.files(tumor_bam_dir, pattern = "\\.bam$", full.names = FALSE) %>%
    str_remove("\\.bam$")

cat(glue("Found {length(tumor_samples)} tumor samples\n"))

# Load sample sheet if using matched normals
if (use_matched_normals) {
    sample_pairs <- read_csv(sample_sheet, col_types = cols(.default = "c")) %>%
        filter(tumor_id %in% tumor_samples)
    
    cat(glue("Found {nrow(sample_pairs)} tumor-normal pairs in sample sheet\n"))
    
    # Validate that normal BAMs exist
    missing_normals <- sample_pairs %>%
        mutate(normal_bam = fs::path(normal_bam_dir, glue("{normal_id}.bam"))) %>%
        filter(!fs::file_exists(normal_bam)) %>%
        pull(normal_id)
    
    if (length(missing_normals) > 0) {
        cat(glue("Warning: Missing normal BAMs for: {paste(missing_normals, collapse = ', ')}\n"))
    }
}

# Function to run snp-pileup with PON
run_snp_pileup_with_pon <- function(tumor_id, pon_file, output_file) {
    tumor_bam <- fs::path(tumor_bam_dir, glue("{tumor_id}.bam"))
    
    if (!fs::file_exists(tumor_bam)) {
        cat(glue("Warning: Tumor BAM not found: {tumor_bam}\n"))
        return(FALSE)
    }
    
    # Read existing PON and append tumor data
    # This approach reads the PON and creates a new pileup with tumor
    cmd <- glue("snp-pileup -g -q15 -Q20 -P100 -r25,0 -d {pon_file} {snp_vcf} {output_file} {tumor_bam}")
    
    cat(glue("Running: {cmd}\n"))
    result <- system(cmd)
    
    return(result == 0 && fs::file_exists(output_file))
}

# Function to run snp-pileup with matched normal
run_snp_pileup_with_normal <- function(tumor_id, normal_id, output_file) {
    tumor_bam <- fs::path(tumor_bam_dir, glue("{tumor_id}.bam"))
    normal_bam <- fs::path(normal_bam_dir, glue("{normal_id}.bam"))
    
    if (!fs::file_exists(tumor_bam)) {
        cat(glue("Warning: Tumor BAM not found: {tumor_bam}\n"))
        return(FALSE)
    }
    
    if (!fs::file_exists(normal_bam)) {
        cat(glue("Warning: Normal BAM not found: {normal_bam}\n"))
        return(FALSE)
    }
    
    cmd <- glue("snp-pileup -g -q15 -Q20 -P100 -r25,0 {snp_vcf} {output_file} {normal_bam} {tumor_bam}")
    
    cat(glue("Running: {cmd}\n"))
    result <- system(cmd)
    
    return(result == 0 && fs::file_exists(output_file))
}

# Function to run FACETS analysis with cnv_facets compatible output
run_facets_analysis <- function(pileup_file, sample_id, output_dir) {
    cat(glue("Processing {sample_id}\n"))
    
    tryCatch({
        # Read pileup data
        rcmat <- readSnpMatrix(pileup_file)
        
        # Pre-process
        xx <- preProcSample(rcmat, gbuild = genome_build)
        
        # Segmentation
        oo <- procSample(xx, cval = cval)
        
        # Fit copy number
        fit <- emcncf(oo)
        
        # Create output files matching cnv_facets format
        sample_output_dir <- fs::path(output_dir, sample_id)
        fs::dir_create(sample_output_dir, recurse = TRUE)
        
        # 1. Save main results (RData format like cnv_facets)
        save(rcmat, xx, oo, fit, file = fs::path(sample_output_dir, glue("{sample_id}.RData")))
        
        # 2. Export segments in cnv_facets format (.seg file)
        if (!is.null(fit$cncf)) {
            segments_df <- fit$cncf %>%
                as_tibble() %>%
                mutate(
                    ID = sample_id,
                    chrom = paste0("chr", chrom),
                    loc.start = start,
                    loc.end = end,
                    num.mark = num.mark,
                    seg.mean = cnlr.median,
                    tcn = tcn.em,
                    lcn = lcn.em,
                    cf = cf.em %||% NA_real_
                ) %>%
                select(ID, chrom, loc.start, loc.end, num.mark, seg.mean, tcn, lcn, cf) %>%
                write_tsv(fs::path(sample_output_dir, glue("{sample_id}.seg")))
            
            # 3. Export copy number calls in cnv_facets format (.cnv file)
            cnv_calls <- fit$cncf %>%
                as_tibble() %>%
                filter(!is.na(tcn.em)) %>%
                mutate(
                    sample = sample_id,
                    chromosome = paste0("chr", chrom),
                    start = start,
                    end = end,
                    copy_number = tcn.em,
                    minor_cn = lcn.em,
                    major_cn = tcn.em - lcn.em,
                    cellular_fraction = cf.em %||% NA_real_,
                    cnv_type = case_when(
                        tcn.em > 2.5 ~ "AMP",
                        tcn.em < 1.5 ~ "DEL", 
                        lcn.em == 0 ~ "LOH",
                        TRUE ~ "NEUTRAL"
                    )
                ) %>%
                select(sample, chromosome, start, end, copy_number, minor_cn, major_cn, 
                       cellular_fraction, cnv_type) %>%
                write_tsv(fs::path(sample_output_dir, glue("{sample_id}.cnv")))
        }
        
        # 4. Generate plots
        pdf(fs::path(sample_output_dir, glue("{sample_id}_spider.pdf")), width = 12, height = 8)
        plotSample(x = oo, emfit = fit)
        dev.off()
        
        # 5. Generate summary file
        summary_stats <- tibble(
            sample = sample_id,
            purity = fit$purity %||% NA_real_,
            ploidy = fit$ploidy %||% NA_real_,
            dipLogR = fit$dipLogR %||% NA_real_,
            loglik = fit$loglik %||% NA_real_,
            flags = if_else(is.null(fit$flags), "", paste(fit$flags, collapse = ",")),
            version = packageVersion("facets") %>% as.character()
        ) %>%
        write_tsv(fs::path(sample_output_dir, glue("{sample_id}_summary.txt")))
        
        return(summary_stats)
        
    }, error = function(e) {
        cat(glue("Error processing {sample_id}: {e$message}\n"))
        return(tibble(sample = sample_id, purity = NA, ploidy = NA, dipLogR = NA, loglik = NA, 
                     flags = "ERROR", version = NA))
    })
}

# Process samples based on approach
if (use_pon) {
    # PON approach - process all tumor samples
    all_summaries <- tumor_samples %>%
        set_names() %>%
        map_dfr(~ {
            cat(glue("\n=== Processing tumor sample: {.x} (PON mode) ===\n"))
            
            # Create pileup file
            pileup_file <- fs::path(output_dir, glue("{.x}_pileup.csv.gz"))
            
            # Run snp-pileup with PON
            success <- run_snp_pileup_with_pon(.x, pon_file, pileup_file)
            
            if (success && fs::file_exists(pileup_file)) {
                run_facets_analysis(pileup_file, .x, output_dir)
            } else {
                cat(glue("Failed to create pileup file for {.x}\n"))
                tibble(sample = .x, purity = NA, ploidy = NA, dipLogR = NA, loglik = NA,
                      flags = "PILEUP_FAILED", version = NA)
            }
        })
        
} else if (use_matched_normals) {
    # Matched normal approach - process tumor-normal pairs
    all_summaries <- sample_pairs %>%
        pmap_dfr(function(tumor_id, normal_id, ...) {
            cat(glue("\n=== Processing tumor-normal pair: {tumor_id} vs {normal_id} ===\n"))
            
            # Create pileup file
            pileup_file <- fs::path(output_dir, glue("{tumor_id}_pileup.csv.gz"))
            
            # Run snp-pileup with matched normal
            success <- run_snp_pileup_with_normal(tumor_id, normal_id, pileup_file)
            
            if (success && fs::file_exists(pileup_file)) {
                run_facets_analysis(pileup_file, tumor_id, output_dir)
            } else {
                cat(glue("Failed to create pileup file for {tumor_id}\n"))
                tibble(sample = tumor_id, purity = NA, ploidy = NA, dipLogR = NA, loglik = NA,
                      flags = "PILEUP_FAILED", version = NA)
            }
        })
}

# Save combined results
if (nrow(all_summaries) > 0) {
    all_summaries %>%
        mutate(
            analysis_type = if_else(use_pon, "PON", "matched_normal"),
            processed_date = Sys.time()
        ) %>%
        rename(sample_id = sample, tumor_purity = purity, genome_ploidy = ploidy) %>%
        write_tsv(fs::path(output_dir, "facets_summary.txt"))
    
    # Create manifest
    manifest <- tibble(
        sample_id = all_summaries$sample,
        analysis_type = if_else(use_pon, "PON", "matched_normal"),
        seg_file = fs::path(output_dir, all_summaries$sample, glue("{all_summaries$sample}.seg")),
        cnv_file = fs::path(output_dir, all_summaries$sample, glue("{all_summaries$sample}.cnv")),
        rdata_file = fs::path(output_dir, all_summaries$sample, glue("{all_summaries$sample}.RData"))
    ) %>%
    write_tsv(fs::path(output_dir, "sample_manifest.txt"))
    
    # Save analysis metadata
    metadata <- tibble(
        analysis_type = if_else(use_pon, "PON", "matched_normal"),
        pon_file = pon_file %||% NA_character_,
        normal_bam_dir = normal_bam_dir %||% NA_character_,
        sample_sheet = sample_sheet %||% NA_character_,
        genome_build = genome_build,
        cval = cval,
        total_samples = nrow(all_summaries),
        successful_samples = sum(!is.na(all_summaries$purity)),
        run_date = Sys.time()
    ) %>%
    write_tsv(fs::path(output_dir, "analysis_metadata.txt"))
}

cat(glue("\nFACETS analysis completed. Results saved in: {output_dir}\n"))
cat(glue("Analysis type: {if_else(use_pon, 'PON', 'Matched Normal')}\n"))
cat(glue("Processed {nrow(all_summaries)} samples\n"))
