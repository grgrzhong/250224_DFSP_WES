#!/bin/bash

#############################################################################
# Initialize the environment
#############################################################################
# Load general configurations
source /lustre1/g/path_my/250224_DFSP_WES/config/env_config.sh
source /lustre1/g/path_my/250224_DFSP_WES/modules/utils/utils.sh

# Log the start time
start_time=$(date +%s)
log_message "Starting variant annotation using GATK Funcotator ..."

# Reference directories/files
export REFERENCE=$REFERENCE_DIR/Gencode/gencode.hg38.v36.primary_assembly.fa
export GERMLINE=$REFERENCE_DIR/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
export ANNOTATION_FILE=$REFERENCE_DIR/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
export INTERVAL=$REFERENCE_DIR/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed

# Input/Output directories/files
export BAM=$BASE_DIR/$PROJECT_DIR/$WORK_DIR/BAM
export PON=$BASE_DIR/$PROJECT_DIR/$WORK_DIR/PON-Mutect/pon.vcf.gz

export OUT_DIR=$BASE_DIR/$PROJECT_DIR/$WORK_DIR/Mutect-Call && mkdir -p $OUT_DIR

# Create logs directory
export LOG_DIR=$OUT_DIR/logs && mkdir -p $LOG_DIR

# Sample sheet file
export SAMPLE_SHEET=$BASE_DIR/$PROJECT_DIR/$WORK_DIR/samplesheet.csv

#############################################################################
# Function to annotate variants using Funcotator
#############################################################################
annotate_with_funcotator() {
    local tumour=$1
    local var_out="$OUT_DIR/${tumour}"
    
    log_message "Annotating variants with Funcotator for sample: ${tumour}"
    
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

    less -S $var_out/${tumour}_annotated.maf.gz | grep -v "#" > $var_out/${tumour}_annotated.tsv
    
    log_message "Completed Funcotator annotation for: ${tumour}"
}

export -f annotate_with_funcotator

# Extract tumor sample IDs
tumor_samples=$(tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f2 | grep -E "\-T")

# Run Funcotator annotation in parallel
echo "$tumor_samples" | parallel -j $PARALLEL_JOBS \
    --joblog "$LOG_DIR/parallel_annotate_funcotator.log" \
    annotate_with_funcotator {}

# Log completion time
log_time_usage $start_time

# Deactivate conda environment if it exists
if [[ -n "${CONDA_DEFAULT_ENV}" ]]; then
    conda deactivate
    log_message "Conda environment deactivated"
fi
