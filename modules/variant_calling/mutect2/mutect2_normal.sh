#!/bin/bash

###############################################################################
# Mutect2 Normal Sample Processing Script
###############################################################################
# Description: This script runs GATK Mutect2 in tumor-only mode on normal samples
#              to generate VCFs for Panel of Normals (PoN) creation. It is the
#              first step in creating a PoN for somatic variant filtering.
#
# Input:  - Normal sample BAM files
# Output: - VCF files for each normal sample
###############################################################################
mutect2_normal() {
    # Validate required global variables
    for var in REFERENCE BAM INTERVAL PON_OUT AVAIL_MEM LOG_DIR; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done
    
    # Parse input parameters
    local normal=$1
    local bam_file="$BAM/$normal/${normal}_recalibrated.bam"
    local output_vcf="$PON_OUT/${normal}_pon.vcf.gz"
    
    # Verify normal BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        echo "ERROR: BAM file not found for $normal: $bam_file" >&2
        return 1
    fi
    
    # Verify reference files exist
    if [[ ! -f "$REFERENCE" ]]; then
        echo "ERROR: Reference genome not found at $REFERENCE" >&2
        return 1
    fi
    
    if [[ ! -f "$INTERVAL" ]]; then
        echo "ERROR: Interval file not found at $INTERVAL" >&2
        return 1
    fi

    # Log start of processing
    log_message "Processing normal sample: $normal"
    
    # Run Mutect2 in tumor-only mode
    gatk --java-options "-Xmx${AVAIL_MEM}g" Mutect2 \
        -R "$REFERENCE" \
        -I "$bam_file" \
        --max-mnp-distance 0 \
        -L "$INTERVAL" \
        --native-pair-hmm-threads 8 \
        -O "$output_vcf" \
        >& "$LOG_DIR/${normal}_Mutect2Normal.log"
    
    # Check if command was successful
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Mutect2 failed for $normal" >&2
        return 1
    fi
    
    # Verify output file exists and is valid
    if [[ ! -f "$output_vcf" ]]; then
        echo "ERROR: Output VCF not created for $normal" >&2
        return 1
    fi
    
    # Verify the output file is a valid VCF
    if ! bcftools view -h "$output_vcf" &>/dev/null; then
        echo "ERROR: Generated VCF for $normal appears to be invalid" >&2
        return 1
    fi
    
    log_message "Successfully processed normal sample: $normal"

}

export -f mutect2_normal
