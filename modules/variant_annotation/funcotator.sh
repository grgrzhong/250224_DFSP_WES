#!/bin/bash

###############################################################################
# Variant Annotation with Funcotator Script
###############################################################################
# Description: This script annotates somatic variants using GATK Funcotator
#
# Input:  - Filtered VCF files from Mutect2
# Output: - Annotated MAF files
#         - Annotated TSV files
###############################################################################

annotate_with_funcotator() {
    # Validate required global variables
    for var in MUTECT_CALL REFERENCE ANNOTATION_FILE INTERVAL; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    local tumour=$1
    local var_out="${MUTECT_CALL}/${tumour}"
    
    # Ensure output directory exists
    mkdir -p "$var_out"
    
    log_message "Annotating variants with Funcotator for sample: ${tumour}"
    
    # Check if input VCF exists
    if [[ ! -f "$var_out/${tumour}_normalized_filtered.vcf.gz" ]]; then
        log_message "ERROR: Input VCF file not found for ${tumour}: $var_out/${tumour}_normalized_filtered.vcf.gz"
        return 1
    fi
    
    # Run Funcotator annotation
    gatk Funcotator \
        -R $REFERENCE \
        -V $var_out/${tumour}_normalized_filtered.vcf.gz \
        -O $var_out/${tumour}_annotated.maf.gz \
        -L $INTERVAL \
        --output-file-format MAF \
        --data-sources-path $ANNOTATION_FILE \
        --ref-version hg38 \
        --remove-filtered-variants true \
        >& ${var_out}/${tumour}.Funcotator.log
    
    # Check if Funcotator ran successfully
    if [[ $? -ne 0 || ! -f "$var_out/${tumour}_annotated.maf.gz" ]]; then
        log_message "ERROR: Funcotator annotation failed for ${tumour}"
        return 1
    fi
    
    # Create a TSV for easier viewing
    log_message "Creating TSV file for ${tumour}"
    less -S $var_out/${tumour}_annotated.maf.gz | grep -v "#" > $var_out/${tumour}_annotated.tsv
    
    log_message "Completed Funcotator annotation for: ${tumour}"
    return 0
}

# Export the function so it's available to GNU parallel
export -f annotate_with_funcotator
