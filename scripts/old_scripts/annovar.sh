#!/bin/bash

###############################################################################
# Variant Annotation with ANNOVAR Script
###############################################################################
# Description: This script annotates somatic variants using ANNOVAR
#
# Input:  - Filtered VCF files from Mutect2
# Output: - Annotated ANNOVAR multianno files
#         - Simplified text files with extracted variant information
###############################################################################

annotate_with_annovar() {
    # Validate required global variables
    for var in  MUTECT_CALL; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    local tumour=$1
    local var_out="${MUTECT_CALL}/${tumour}"
    local var_ann=${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/variant_annotation/annovar/${tumour}
    # Ensure output directory exists
    mkdir -p "$var_out"
    mkdir -p "$var_ann"

    log_message "Processing ANNOVAR annotation for sample: ${tumour}"
    
    # Check if input VCF exists
    if [[ ! -f "$var_out/${tumour}_normalized_filtered.vcf.gz" ]]; then
        log_message "ERROR: Input VCF file not found for ${tumour}: $var_out/${tumour}_normalized_filtered.vcf.gz"
        return 1
    fi
    
    # Run ANNOVAR annotation
    perl $BASE_DIR/Software/annovar//table_annovar.pl $var_out/${tumour}_normalized_filtered.vcf.gz \
        $BASE_DIR/Software/annovar/humandb/ \
        -buildver hg38 -out $var_ann/${tumour} -remove \
        -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
        -operation gx,r,f,f,f,f,f \
        -nastring . -polish -xreffile $BASE_DIR/Software/annovar/example/gene_fullxref.txt \
        --otherinfo --vcfinput \
        >& ${var_ann}/${tumour}.Annovar.log
    
    # Check if ANNOVAR ran successfully
    if [[ $? -ne 0 || ! -f "$var_ann/${tumour}.hg38_multianno.txt" ]]; then
        log_message "ERROR: ANNOVAR annotation failed for ${tumour}"
        return 1
    fi
    
    # Extract additional fields from the VCF and create a simplified output
    log_message "Creating simplified annotation file for ${tumour}"
    less -S $var_ann/${tumour}.hg38_multianno.txt | \
        awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' \
        > $var_ann/${tumour}_annovar.txt
    
    # Check if output file was created successfully
    if [[ ! -f "$var_ann/${tumour}_annovar.txt" ]]; then
        log_message "ERROR: Failed to create simplified annotation file for ${tumour}"
        return 1
    fi
    
    log_message "Completed ANNOVAR annotation for: ${tumour}"
    return 0
}

# Export the function so it's available to GNU parallel
export -f annotate_with_annovar
