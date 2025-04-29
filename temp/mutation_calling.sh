#!/bin/bash

#############################################################################
# Workflow: Somatic Variant Discovery Pipeline
# Description: Orchestrates variant calling modules for WES analysis
#############################################################################

# Load the project configuration
source "/lustre1/g/path_my/250224_DFSP_WES/conf/project_config.sh"

# Print project settings
echo "====================== Project Settings ======================"
echo "BASE_DIR:              ${BASE_DIR}"
echo "PROJECT_DIR:           ${PROJECT_DIR}"
echo "WORK_DIR:              ${WORK_DIR}"
echo "REFERENCE_DIR:         ${REFERENCE_DIR}"
echo "SAMPLE_SHEET:          ${SAMPLE_SHEET}"
echo "Total cases:           ${CASE_COUNT}"
echo "Total Tumour samples:  ${TUMOUR_COUNT}"
echo "Total Normal samples:  ${NORMAL_COUNT}"
echo "Parallel jobs:         ${PARALLEL_JOBS}"
echo "==============================================================="

#############################################################################
# Load conda environment (modify this path according to your HPC setup)
#############################################################################
# if [ -f "/home/zhonggr/miniforge3/etc/profile.d/conda.sh" ]; then
#     source "/home/zhonggr/miniforge3/etc/profile.d/conda.sh"
# else
#     echo "Error: Conda initialization file not found"
#     exit 1
# fi

# conda activate varcall || echo "Error: Failed to activate conda environment 'varcall'"

#############################################################################
# Mutect2 call
#############################################################################
run_Mutect2Normal() {
    # Define local variables for paths and sample data
    local $NORMAL_SAMPLES="${NORMAL_SAMPLES}"
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
    export BAM="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/BAM"
    export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
    export PON_OUT="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/PON-Mutect"
    export AVAIL_MEM="${AVAIL_MEM}"

    # Create logs directory if not exists
    LOG_DIR="${PON_OUT}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/mutect2_normal.sh"

    # Run GATK Mutect2 in parallel
    echo "$NORMAL_SAMPLES" | parallel \
        --env REFERENCE,INTERVAL,BAM,PON_OUT,AVAIL_MEM \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_Mutect2Normal.log" \
        mutect2_normal {}
}

run_GenomicsDBImport() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export BAM="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/BAM"
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
    export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
    export AVAIL_MEM="${AVAIL_MEM}"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/genomics_db_import.sh"

    # Run GATK GenomicsDBImport in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env BAM,REFERENCE,INTERVAL,AVAIL_MEM \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_GenomicsDBImport.log" \
        genomics_db_import {}
}

run_CreateSomaticPanelOfNormals() {
    # Define local variables for paths and sample data
    export PON_OUT="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/PON-Mutect"
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
    export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
    export AVAIL_MEM="${AVAIL_MEM}"

    # Create logs directory if not exists
    LOG_DIR="${PON_OUT}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/create_pon.sh"

    # Run GATK CreateSomaticPanelOfNormals in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env PON_OUT,REFERENCE,INTERVAL.AVAIL_MEM \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_CreateSomaticPanelOfNormals.log" \
        create_pon {}
}

run_GetpileupSummaries() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"
    export BAM="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/BAM"
    export REFERENCE_DIR="${REFERENCE_DIR}"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/get_pileup_summaries.sh"

    # Run GATK GetpileupSummaries in parallel
    echo "$TUMOUR_SAMPLES" | parallel --env BAM,MUTECT_CALL,REFERENCE_DIR \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_GetpileupSummaries.log" \
        get_pileup_summaries {}
}

run_CalculateContamination() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/calculate_contamination.sh"

    # Run GATK CalculateContamination in parallel
    echo "$TUMOUR_SAMPLES" | parallel --env MUTECT_CALL \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_CalculateContamination.log" \
        calculate_contamination {}
}

run_Mutect2CallVariant() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export BAM="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/BAM"
    export PON="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/PON-Mutect/pon.vcf.gz"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"
    
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
    export GERMLINE="${REFERENCE_DIR}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
    export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
    export AVAIL_MEM="${AVAIL_MEM}"
    
    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/parallel_logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/mutect2_call.sh"

    # Run Mutect2 in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env BAM,PON,MUTECT_CALL,REFERENCE,GERMLINE,INTERVAL \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_Mutect2CallVariant.log" \
        mutect2_call {}
}

run_LearnReadOrientationModel() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/learn_read_orientation_model.sh"

    # Run GATK LearnReadOrientationModel in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env MUTECT_CALL \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_LearnReadOrientationModel.log" \
        learn_read_orientation_model {}
}

run_FilterMutectCalls() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
    export AVAIL_MEM="${AVAIL_MEM}"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/filter_mutect_calls.sh"

    # Run GATK FilterMutectCalls in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env MUTECT_CALL,REFERENCE \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_FilterMutectCalls.log" \
        filter_mutect_calls {}
}

run_NormalizeReads() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/normalize_reads.sh"

    # Run GATK NormalizeReads in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env MUTECT_CALL,REFERENCE \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_NormalizeReads.log" \
        normalize_reads {}
}

run_FuncotatorAnnotation() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"
    export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
    export ANNOTATION_FILE="${REFERENCE_DIR}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/"
    export INTERVAL="${REFERENCE_DIR}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/funcotator.sh"

    # Run GATK Funcotator annotation in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env MUTECT_CALL,REFERENCE,ANNOTATION_FILE,INTERVAL \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_Funcotator.log" \
        annotate_with_funcotator {}
}

run_AnnovarAnnotation() {
    # Define local variables for paths and sample data
    local TUMOUR_SAMPLES="${TUMOUR_SAMPLES}"
    export MUTECT_CALL="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"

    # Create logs directory if not exists
    LOG_DIR="${MUTECT_CALL}/logs" && mkdir -p "${LOG_DIR}"

    # Load module
    source "${BASE_DIR}/${PROJECT_DIR}/scripts/annovar.sh"

    # Run ANNOVAR annotation in parallel
    echo "$TUMOUR_SAMPLES" | parallel \
        --env MUTECT_CALL \
        --jobs $PARALLEL_JOBS \
        --joblog "$LOG_DIR/parallel_Annovar.log" \
        annotate_with_annovar {}
}

# # Main Workflow
# main() {
#     # Validate environment
#     check_environment || exit 1
    
#     # Read sample sheet
#     local sample_sheet="${DATA_DIR}/study_samples.csv"
#     [[ ! -f "$sample_sheet" ]] && {
#         log_error "Sample sheet not found: $sample_sheet"
#         exit 1
#     }
    
#     # Process each sample pair
#     while IFS=, read -r patient tumor normal; do
#         [[ "$patient" == "patient" ]] && continue  # Skip header
        
#         # Create output directories
#         mkdir -p "${OUT_DIR}/${tumor}"
        
#         # Run workflow steps
#         run_preprocessing "$tumor"
#         run_preprocessing "$normal"
#         run_alignment "$tumor"
#         run_alignment "$normal"
#         run_variant_calling "$tumor" "$normal"
#     done < "$sample_sheet"
    
#     log_message "Workflow completed successfully"
# }

# # Execute workflow
# main "$@"