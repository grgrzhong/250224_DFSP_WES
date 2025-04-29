#!/bin/bash

###############################################################################
# Call Somatic Variants Script
###############################################################################
# Description: This script performs somatic variant calling using GATK Mutect2
#              for tumor-normal or tumor-only samples.
#
# Input:  - Recalibrated BAM files
# Output: - Raw VCF file with somatic variant calls
#         - Realigned BAM file
###############################################################################

mutect2_call() {
    # Validate required global variables
    for var in BAM PON MUTECT_CALL REFERENCE INTERVAL GERMLINE; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    local tumour=$1
    local case=$(echo "$tumour" | sed -E 's/^(.*?)-(T|N).*$/\1/')
    local normal="${case}-N"
    local out_dir="$MUTECT_CALL/${tumour}"
    
    local avail_mem=${AVAIL_MEM:-4}
    # Ensure output directory exists
    mkdir -p "$out_dir"
    
    log_message "Starting variants calling using Mutect2 ..."

    # Check for normal sample
    if [[ -d "${BAM}/${normal}" && -f "${BAM}/${normal}/${normal}_recalibrated.bam" ]]; then
        
        log_message "Processing ${tumour} with matched normal sample ${normal} ..."
        
        gatk --java-options "-Xmx${avail_mem}g" Mutect2 \
            -I "${BAM}/${tumour}/${tumour}_recalibrated.bam" \
            -I "${BAM}/${normal}/${normal}_recalibrated.bam" \
            -normal "$normal" \
            -R "$REFERENCE" \
            -L "$INTERVAL" \
            --germline-resource "$GERMLINE" \
            --panel-of-normals "$PON" \
            --f1r2-tar-gz "${out_dir}/${tumour}.f1r2.tar.gz" \
            --native-pair-hmm-threads 8 \
            --callable-depth 20 \
            -O "${out_dir}/${tumour}_unfiltered.vcf.gz" \
            -bamout "${out_dir}/${tumour}_realigned.bam" \
            >& "${out_dir}/${tumour}.Mutect2Call.log"
    else
        
        log_message "Processing ${tumour} without matched normal sample ..."
        
        gatk --java-options "-Xmx${avail_mem}g" Mutect2 \
            -I "${BAM}/${tumour}/${tumour}_recalibrated.bam" \
            -R "$REFERENCE" \
            -L "$INTERVAL" \
            --germline-resource "$GERMLINE" \
            --panel-of-normals "$PON" \
            --f1r2-tar-gz "${out_dir}/${tumour}.f1r2.tar.gz" \
            --native-pair-hmm-threads 8 \
            --callable-depth 20 \
            -O "${out_dir}/${tumour}_unfiltered.vcf.gz" \
            -bamout "${out_dir}/${tumour}_realigned.bam" \
            >& "${out_dir}/${tumour}.Mutect2Call.log"
    fi
}

export -f mutect2_call