#!/bin/bash

###############################################################################
# GenomicsDB Import Module
###############################################################################
# Description: This module imports normal sample VCFs into a GenomicsDB workspace
#              for Panel of Normals (PoN) creation. It handles memory allocation,
#              command construction, and execution for GATK GenomicsDBImport.
#
# Input:  - Normal sample VCFs (.vcf.gz)
# Output: - GenomicsDB workspace for Panel of Normals
###############################################################################

import_genomics_db() {
    # Validate required global variables
    for var in REFERENCE INTERVAL PON_OUT; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    # Verify reference files exist
    if [[ ! -f "$REFERENCE" ]]; then
        echo "ERROR: Reference genome not found at $REFERENCE" >&2
        return 1
    fi
    
    if [[ ! -f "$INTERVAL" ]]; then
        echo "ERROR: Interval file not found at $INTERVAL" >&2
        return 1
    fi

    # Create GenomicsDB workspace
    cmd="gatk --java-options '-Xmx${AVAIL_MEM}g' GenomicsDBImport \
        -R $REFERENCE \
        -L $INTERVAL \
        --genomicsdb-workspace-path $PON_OUT/pon_db \
        --merge-input-intervals true"

    log_message "Importing indexed VCF files into GenomicsDB"

    for vcf_file in $(ls $PON_OUT/*.vcf.gz); do

        cmd+=" -V $vcf_file"

    done

    eval $cmd

    # Verify the workspace was created
    if [[ ! -d "$PON_OUT/pon_db" ]]; then
        echo "ERROR: GenomicsDB workspace not created" >&2
        return 1
    fi

    log_message "Successfully created GenomicsDB workspace at $PON_OUT/pon_db"

    return 0
}

export -f import_genomics_db