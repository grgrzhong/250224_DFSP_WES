#!/bin/bash

###############################################################################
# Normalize Reads Script
###############################################################################
# Description: This script normalizes variants in a VCF file using bcftools norm
#              to split multiallelic sites and left-align indels. It also filters
#              for PASS variants and creates a tabix index.
#
# Input:  - Filtered VCF from Mutect2
# Output: - Normalized and filtered VCF file with index
###############################################################################

normalize_reads() {
    # Validate required global variables
    for var in MUTECT_CALL REFERENCE; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    # Parse input parameters
    local tumour=$1
    local out_dir="$MUTECT_CALL/${tumour}"
    local input_vcf="${out_dir}/${tumour}_filtered.vcf.gz"
    local temp_vcf="${out_dir}/${tumour}_normalized.vcf.gz"
    local output_vcf="${out_dir}/${tumour}_normalized_filtered.vcf.gz"
    
    # Check if output directory exists
    if [[ ! -d "$out_dir" ]]; then
        echo "ERROR: Output directory not found: $out_dir" >&2
        return 1
    fi
    
    # Check if input VCF exists
    if [[ ! -f "$input_vcf" ]]; then
        echo "ERROR: Input VCF not found: $input_vcf" >&2
        return 1
    fi
    
    # Check if reference genome exists
    if [[ ! -f "$REFERENCE" ]]; then
        echo "ERROR: Reference genome not found at $REFERENCE" >&2
        return 1
    fi
    
    log_message "Normalizing variants for sample: ${tumour}"
    
    # Step 1: Normalize variants (split multi-allelics and left-align indels)
    bcftools norm -m-both -f "$REFERENCE" -Oz \
        -o "$temp_vcf" "$input_vcf"
        
    # Check if normalization was successful
    if [[ $? -ne 0 || ! -f "$temp_vcf" ]]; then
        echo "ERROR: bcftools norm failed for ${tumour}" >&2
        return 1
    fi
    
    # Step 2: Filter for PASS variants
    bcftools view -f PASS "$temp_vcf" \
        -o "$output_vcf"
        
    # Check if filtering was successful
    if [[ $? -ne 0 || ! -f "$output_vcf" ]]; then
        echo "ERROR: bcftools view failed for ${tumour}" >&2
        return 1
    fi
    
    # Step 3: Create index using tabix
    tabix "$output_vcf"
    
    # Check if indexing was successful
    if [[ $? -ne 0 || ! -f "${output_vcf}.tbi" ]]; then
        echo "ERROR: tabix indexing failed for ${tumour}" >&2
        return 1
    fi
    
    # Step 4: Clean up temporary files
    if [[ -f "$temp_vcf" ]]; then
        rm -f "$temp_vcf"
        log_message "Removing temporary file: $(basename "$temp_vcf")"
    fi
    
    log_message "Completed normalization for: ${tumour}"

}

export -f normalize_reads