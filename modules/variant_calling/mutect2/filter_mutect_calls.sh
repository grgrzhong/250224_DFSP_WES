#!/bin/bash

###############################################################################
# Filter Mutect Calls Script
###############################################################################
# Description: This script filters somatic variant calls from Mutect2 using GATK
#              FilterMutectCalls. It applies orientation bias filtering, 
#              contamination filtering, and other quality metrics.
#
# Input:  - Unfiltered VCF from Mutect2
#         - Read orientation model
#         - Contamination table
#         - Segmentation table
# Output: - Filtered VCF file
###############################################################################

filter_mutect_calls() {
    # Validate required global variables
    for var in MUTECT_CALL REFERENCE AVAIL_MEM; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    local tumour=$1
    local out_dir="$MUTECT_CALL/${tumour}"
    local case=$(echo "$tumour" | sed -E 's/^(.*?)-(T|N).*$/\1/')
    
    # Check if output directory exists
    if [[ ! -d "$out_dir" ]]; then
        echo "ERROR: Output directory not found: $out_dir" >&2
        return 1
    fi
    
    # Check if required input files exist
    if [[ ! -f "${out_dir}/${tumour}_unfiltered.vcf.gz" ]]; then
        echo "ERROR: Unfiltered VCF not found for ${tumour}" >&2
        return 1
    fi
    
    if [[ ! -f "${out_dir}/${tumour}.read-orientation-model.tar.gz" ]]; then
        echo "ERROR: Read orientation model not found for ${tumour}" >&2
        return 1
    fi
    
    if [[ ! -f "${out_dir}/${tumour}.contamination.table" ]]; then
        echo "ERROR: Contamination table not found for ${tumour}" >&2
        return 1
    fi
    
    if [[ ! -f "${out_dir}/${tumour}.segments.table" ]]; then
        echo "ERROR: Segments table not found for ${tumour}" >&2
        return 1
    fi
    
    # Check if reference genome exists
    if [[ ! -f "$REFERENCE" ]]; then
        echo "ERROR: Reference genome not found at $REFERENCE" >&2
        return 1
    fi
    
    log_message "Starting FilterMutectCalls for sample: ${tumour}"
    
    gatk --java-options "-Xmx${AVAIL_MEM}g" FilterMutectCalls \
        -V "${out_dir}/${tumour}_unfiltered.vcf.gz" \
        -R "$REFERENCE" \
        --ob-priors "${out_dir}/${tumour}.read-orientation-model.tar.gz" \
        --contamination-table "${out_dir}/${tumour}.contamination.table" \
        --tumor-segmentation "${out_dir}/${tumour}.segments.table" \
        --min-allele-fraction 0.01 \
        --unique-alt-read-count 1 \
        --stats "${out_dir}/${tumour}_unfiltered.vcf.gz.stats" \
        -O "${out_dir}/${tumour}_filtered.vcf.gz" \
        >& "${out_dir}/${tumour}.FilterMutectCalls.log"
    
    # Check if the command was successful
    if [[ $? -ne 0 ]]; then
        echo "ERROR: FilterMutectCalls failed for ${tumour}" >&2
        return 1
    fi

    # Verify output file exists
    if [[ ! -f "${out_dir}/${tumour}_filtered.vcf.gz" ]]; then
        echo "ERROR: Filtered VCF not created for ${tumour}" >&2
        return 1
    fi

    log_message "Completed FilterMutectCalls for: ${tumour}"
}

export -f filter_mutect_calls