#!/bin/bash

# Load general configurations
source /lustre1/g/path_my/250224_DFSP_WES/config/env_config.sh
source /lustre1/g/path_my/250224_DFSP_WES/modules/utils.sh

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

# Function to annotate variants using annovar
annotate_with_annovar() {
    local tumour=$1
    local var_out="$OUT_DIR/${tumour}"
    
    log_message "Processing Annovar annotation for sample: ${tumour}"
    
    perl $BASE_DIR/Software/annovar/table_annovar.pl $var_out/${tumour}_normalized_filtered.vcf.gz \
        $BASE_DIR/Software/annovar/humandb/ \
        -buildver hg38 -out $var_out/${tumour} -remove \
        -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
        -operation gx,r,f,f,f,f,f \
        -nastring . -polish -xreffile $BASE_DIR/Software/annovar/example/gene_fullxref.txt \
        --otherinfo --vcfinput \
        >& ${var_out}/${tumour}.Annovar.log

    less -S $var_out/${tumour}.hg38_multianno.txt | \
        awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' \
        > $var_out/${tumour}_annovar.txt
    
    log_message "Completed Annovar annotation for: ${tumour}"
}

export -f annotate_with_annovar

# Extract tumor sample IDs
tumor_samples=$(tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f2 | grep -E "\-T")

# Run Funcotator annotation in parallel
echo "$tumor_samples" | parallel -j $PARALLEL_JOBS \
    --joblog "$LOG_DIR/parallel_annotate_annovar.log" \
    annotate_with_annovar {}

# Log completion time
log_time_usage $start_time

# Deactivate conda environment if it exists
if [[ -n "${CONDA_DEFAULT_ENV}" ]]; then
    conda deactivate
    log_message "Conda environment deactivated"
fi