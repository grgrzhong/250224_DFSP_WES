library(facets)
library(pctGCdata)
library(tidyverse)
library(parallel)
library(glue)

# Configuration
tumor_bam_dir <- "/path/to/tumor/bams" # Update with your tumor BAM directory
normal_bam_dir <- "/path/to/normal/bams" # Update with your normal BAM directory
output_dir <- "/home/zhonggr/projects/250224_DFSP_WES/facets_output"
reference_genome <- "/path/to/reference.fa" # Update with your reference genome
snp_vcf <- "/path/to/dbsnp.vcf" # Update with SNP positions file

# Option 1: Use existing Mutect2 PON
mutect2_pon <- "/path/to/mutect2_pon.vcf.gz" # Your Mutect2 PON file
use_mutect2_pon <- TRUE # Set to FALSE to use BAM-based approach

# Option 2: Use subset of normal BAMs (recommended if Mutect2 PON doesn't work well)
max_normals_for_pon <- 30 # Limit number of normals to avoid memory issues

# Create output directory
fs::dir_create(output_dir, recurse = TRUE)

# Get sample lists using tidyverse approach
tumor_samples <- list.files(tumor_bam_dir, pattern = "\\.bam$", full.names = FALSE) %>%
    str_remove("\\.bam$")

normal_samples <- list.files(normal_bam_dir, pattern = "\\.bam$", full.names = FALSE) %>%
    str_remove("\\.bam$")

cat(glue("Found {length(tumor_samples)} tumor samples and {length(normal_samples)} normal samples\n"))

# Modified function to handle different PON approaches
run_snp_pileup <- function(tumor_id, normal_source, output_file, use_pon = FALSE) {
    # snp-pileup PURPOSE:
    # 1. Extracts read counts for reference (REF) and alternate (ALT) alleles at SNP positions
    # 2. Creates a matrix with tumor vs normal allele frequencies
    # 3. This data is used by FACETS to:
    #    - Calculate Log R ratios (copy number changes)
    #    - Calculate B-allele frequencies (allelic imbalance)
    #    - Detect loss of heterozygosity (LOH)
    #    - Estimate tumor purity and ploidy
    
    tumor_bam <- fs::path(tumor_bam_dir, glue("{tumor_id}.bam"))

    if (use_pon && fs::file_exists(mutect2_pon)) {
        # Option 1: Use Mutect2 PON sites - these are high-confidence germline variants
        # that provide good heterozygous SNPs for allelic balance analysis
        temp_snp_file <- tempfile(fileext = ".vcf")

        # Create proper VCF header and extract SNP positions
        cat("##fileformat=VCFv4.2\n", file = temp_snp_file)
        cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = temp_snp_file, append = TRUE)

        # Extract positions from PON (these are recurrent germline variants)
        cmd_extract <- glue("zcat {mutect2_pon} | grep -v '^#' | awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,\".\",$4,$5,\".\",\"PASS\",\".\"}}' >> {temp_snp_file}")
        system(cmd_extract)

        # Use subset of normal BAMs to create baseline allele frequencies
        selected_normals <- head(normal_source, max_normals_for_pon)
        normal_bam_string <- str_c(selected_normals, collapse = " ")

        # snp-pileup command explanation:
        # -g: discard reads with poor mapping quality
        # -q15: minimum base quality 15
        # -Q20: minimum mapping quality 20  
        # -P100: maximum coverage per position
        # -r25,0: minimum coverage 25, no maximum deletions
        cmd <- glue("snp-pileup -g -q15 -Q20 -P100 -r25,0 {temp_snp_file} {output_file} {normal_bam_string} {tumor_bam}")
        
        cat(glue("Using Mutect2 PON sites with {length(selected_normals)} normal samples\n"))
        
    } else {
        # Option 2: Use dbSNP - comprehensive SNP database for genome-wide coverage
        selected_normals <- head(normal_source, max_normals_for_pon)
        normal_bam_string <- str_c(selected_normals, collapse = " ")

        cmd <- glue("snp-pileup -g -q15 -Q20 -P100 -r25,0 {snp_vcf} {output_file} {normal_bam_string} {tumor_bam}")
        
        cat(glue("Using dbSNP with {length(selected_normals)} normal samples from BAMs\n"))
    }

    # Output format: CSV with columns for each sample showing:
    # Chromosome, Position, REF_count_normal, ALT_count_normal, REF_count_tumor, ALT_count_tumor
    cat(glue("Running: {cmd}\n"))
    result <- system(cmd)

    # Clean up temporary file if created
    if (exists("temp_snp_file") && file.exists(temp_snp_file)) {
        unlink(temp_snp_file)
    }

    return(result == 0 && file.exists(output_file))
}

# Function to create a custom PON from normal BAMs (alternative approach)
create_facets_pon <- function(normal_bams, pon_output_file) {
    cat("Creating FACETS-compatible PON from", length(normal_bams), "normal samples\n")

    # Create a pooled normal by randomly selecting a subset
    selected_normals <- sample(normal_bams, min(20, length(normal_bams)))

    # Create a "fake" tumor-normal pair using two different normals
    if (length(selected_normals) >= 2) {
        fake_tumor <- selected_normals[1]
        fake_normals <- selected_normals[-1]

        normal_bam_string <- paste(fake_normals, collapse = " ")

        cmd <- sprintf(
            "snp-pileup -g -q15 -Q20 -P100 -r25,0 %s %s %s %s",
            snp_vcf, pon_output_file, normal_bam_string, fake_tumor
        )

        system(cmd)
        return(file.exists(pon_output_file))
    }

    return(FALSE)
}

# Function to run FACETS analysis with cnv_facets compatible output
run_facets_analysis <- function(pileup_file, sample_id, output_dir) {
    cat(glue("Processing {sample_id}\n"))

    # Read pileup data
    rcmat <- readSnpMatrix(pileup_file)

    # Pre-process
    xx <- preProcSample(rcmat, gbuild = "hg38") # Change to hg19 if needed

    # Segmentation
    oo <- procSample(xx, cval = 150)

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
                chrom = paste0("chr", chrom),  # Add chr prefix like cnv_facets
                loc.start = start,
                loc.end = end,
                num.mark = num.mark,
                seg.mean = cnlr.median,
                # Additional cnv_facets specific columns
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
                # Classify amplifications/deletions like cnv_facets
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

    # 4. Generate plots in cnv_facets style
    pdf(fs::path(sample_output_dir, glue("{sample_id}_spider.pdf")), width = 12, height = 8)
    plotSample(x = oo, emfit = fit)
    dev.off()

    # 5. Generate summary file matching cnv_facets format
    summary_stats <- tibble(
        sample = sample_id,
        purity = fit$purity %||% NA_real_,
        ploidy = fit$ploidy %||% NA_real_,
        dipLogR = fit$dipLogR %||% NA_real_,
        loglik = fit$loglik %||% NA_real_,
        # Additional metrics like cnv_facets
        flags = if_else(is.null(fit$flags), "", paste(fit$flags, collapse = ",")),
        version = packageVersion("facets") %>% as.character()
    ) %>%
    write_tsv(fs::path(sample_output_dir, glue("{sample_id}_summary.txt")))

    return(summary_stats)
}

# Main workflow
cat("Starting FACETS analysis workflow\n")

# Prepare normal BAM file paths using fs::path
normal_bam_paths <- fs::path(normal_bam_dir, glue("{normal_samples}.bam"))

# Option 3: Create a single PON file for all samples (memory efficient)
if (use_mutect2_pon || length(normal_bam_paths) > max_normals_for_pon) {
    cat("Creating pooled normal reference\n")
    pon_pileup_file <- file.path(output_dir, "pooled_normal_pon.csv.gz")

    if (!file.exists(pon_pileup_file)) {
        create_facets_pon(normal_bam_paths, pon_pileup_file)
    }
}

# Process each tumor sample using map
all_summaries <- tumor_samples %>%
    set_names() %>%
    map_dfr(~ {
        cat(glue("\n=== Processing tumor sample: {.x} ===\n"))

        # Create pileup file
        pileup_file <- fs::path(output_dir, glue("{.x}_pileup.csv.gz"))

        # Run snp-pileup with chosen approach
        success <- run_snp_pileup(.x, normal_bam_paths, pileup_file, use_mutect2_pon)

        if (success && fs::file_exists(pileup_file)) {
            # Run FACETS analysis
            tryCatch(
                {
                    run_facets_analysis(pileup_file, .x, output_dir)
                },
                error = function(e) {
                    cat(glue("Error processing {.x}: {e$message}\n"))
                    tibble(sample_id = .x, purity = NA, ploidy = NA, dipLogR = NA, loglik = NA)
                }
            )
        } else {
            cat(glue("Failed to create pileup file for {.x}\n"))
            tibble(sample_id = .x, purity = NA, ploidy = NA, dipLogR = NA, loglik = NA)
        }
    })

# Save combined summary in cnv_facets format
if (nrow(all_summaries) > 0) {
    # Create master summary file like cnv_facets
    all_summaries %>%
        rename(
            sample_id = sample,
            tumor_purity = purity,
            genome_ploidy = ploidy
        ) %>%
        write_tsv(fs::path(output_dir, "facets_summary.txt"))
    
    # Create a manifest file listing all processed samples
    manifest <- tibble(
        sample_id = all_summaries$sample,
        seg_file = fs::path(output_dir, all_summaries$sample, glue("{all_summaries$sample}.seg")),
        cnv_file = fs::path(output_dir, all_summaries$sample, glue("{all_summaries$sample}.cnv")),
        rdata_file = fs::path(output_dir, all_summaries$sample, glue("{all_summaries$sample}.RData"))
    ) %>%
    write_tsv(fs::path(output_dir, "sample_manifest.txt"))
}

cat(glue("\nFACETS analysis completed. Results saved in: {output_dir}\n"))
