#!/bin/bash

###############################################################################
# Get Pileup Summaries Script
###############################################################################
# Description: This script generates pileup summaries for tumor and matched 
#              normal samples using GATK GetPileupSummaries tool. The pileup 
#              summaries are used forcontamination estimation in variant 
#              calling pipeline.
#
# Input:  - BAM files (recalibrated)
#         - Common variants VCF for contamination estimation
# Output: - Pileup summary tables for tumor and normal samples
###############################################################################

get_pileup_summaries() {
    # Validate required global variables
    for var in BAM MUTECT_CALL REFERENCE_DIR; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    local tumour=$1
    local case=$(echo "$tumour" | sed -E 's/^(.*?)-(T|N).*$/\1/')
    local normal="${case}-N"
    local out_dir="$MUTECT_CALL/${tumour}"

    mkdir -p "$out_dir"
    
    log_message "Starting GetPileupSummaries for ${tumour}"
    
    # Get Pileup Summaries for tumor
    gatk GetPileupSummaries \
        -I $BAM/$tumour/${tumour}_recalibrated.bam \
        -V $REFERENCE_DIR/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -L $REFERENCE_DIR/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -O ${out_dir}/${tumour}.getpileupsummaries.table \
        >& ${out_dir}/${tumour}.GetPileupSummaries.log

    # Get Pileup Summaries for normal if exists
    if [ -d ${BAM}/${normal} ]; then
        
        log_message "Found matched normal sample: ${normal}"
        log_message "Starting GetPileupSummaries for ${normal}"

        gatk GetPileupSummaries \
            -I $BAM/$normal/${normal}_recalibrated.bam \
            -V $REFERENCE_DIR/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -L $REFERENCE_DIR/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -O ${out_dir}/${normal}.getpileupsummaries.table \
            >& ${out_dir}/${normal}.GetPileupSummaries.log
    fi
}

export -f get_pileup_summaries
