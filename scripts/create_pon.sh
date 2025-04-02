#!/bin/bash

###############################################################################
# Panel of Normals Creation Script
###############################################################################
# Description: This script creates a Panel of Normals (PoN) for Mutect2 somatic
#              variant calling. It processes normal samples and generates a VCF
#              containing sites of common technical artifacts or germline variation.
#
# Input:  - Normal sample VCFs or GenomicsDB workspace
# Output: - Panel of Normals VCF file
###############################################################################

create_pon() {
    # Validate required global variables
    for var in REFERENCE INTERVAL PON_OUT AVAIL_MEM LOG_DIR; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done
    
    # Verify required files exist
    if [[ ! -f "$REFERENCE" ]]; then
        echo "ERROR: Reference genome not found at $REFERENCE" >&2
        return 1
    fi
    
    if [[ ! -f "$INTERVAL" ]]; then
        echo "ERROR: Interval file not found at $INTERVAL" >&2
        return 1
    fi
    
    # Verify genomics database exists
    if [[ ! -d "$PON_OUT/pon_db" ]]; then
        echo "ERROR: GenomicsDB not found at $PON_OUT/pon_db" >&2
        echo "       Please run GenomicsDBImport first before creating PoN" >&2
        return 1
    fi
    
    log_message "Starting creation of Panel of Normals..."
    
    # Create Panel of Normals
    gatk --java-options "-Xmx${MAX_MEM}g" CreateSomaticPanelOfNormals \
        -R "$REFERENCE" \
        -L "$INTERVAL" \
        -V "gendb://$PON_OUT/pon_db" \
        -O "$PON_OUT/pon.vcf.gz" \
        >& "$LOG_DIR/CreateSomaticPanelOfNormals.log"
        
    # Check if command was successful
    if [[ $? -ne 0 ]]; then
        echo "ERROR: CreateSomaticPanelOfNormals failed" >&2
        return 1
    fi
    
    # Verify output file exists and is valid
    if [[ ! -f "$PON_OUT/pon.vcf.gz" ]]; then
        echo "ERROR: Output PoN file not created" >&2
        return 1
    fi
    
    # Verify the output file is a valid VCF
    if ! bcftools view -h "$PON_OUT/pon.vcf.gz" &>/dev/null; then
        echo "ERROR: Generated PoN file appears to be invalid" >&2
        return 1
    fi
    
    log_message "Successfully created Panel of Normals at $PON_OUT/pon.vcf.gz"

}

export -f create_pon
