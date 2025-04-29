#!/bin/bash

###############################################################################
# Learn Read Orientation Model Script
###############################################################################
# Description: This script creates a read orientation model using GATK
#              LearnReadOrientationModel tool. The model is used to filter
#              orientation bias artifacts during Mutect2 filtering.
#
# Input:  - F1R2 counts from Mutect2
# Output: - Read orientation model file (.tar.gz)
###############################################################################

learn_read_orientation_model() {
    # Validate required global variables
    for var in MUTECT_CALL; do
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
    
    # Check if F1R2 tar.gz file exists
    if [[ ! -f "${out_dir}/${tumour}.f1r2.tar.gz" ]]; then
        echo "ERROR: F1R2 file not found for ${tumour}" >&2
        return 1
    fi

    log_message "Learning Read Orientation Model for sample: ${tumour}"
    
    gatk LearnReadOrientationModel \
        -I "${out_dir}/${tumour}.f1r2.tar.gz" \
        -O "${out_dir}/${tumour}.read-orientation-model.tar.gz" \
        >& "${out_dir}/${tumour}.LearnReadOrientationModel.log"
        
    # Check if the command was successful
    if [[ $? -ne 0 ]]; then
        echo "ERROR: LearnReadOrientationModel failed for ${tumour}" >&2
        return 1
    fi
    
    # Verify output file exists
    if [[ ! -f "${out_dir}/${tumour}.read-orientation-model.tar.gz" ]]; then
        echo "ERROR: Output orientation model not created for ${tumour}" >&2
        return 1
    fi
    
    log_message "Completed Read Orientation Model for: ${tumour}"
}

export -f learn_read_orientation_model