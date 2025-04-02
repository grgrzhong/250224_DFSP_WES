#!/bin/bash

###############################################################################
# Strelka2 Somatic Variant Calling Module
###############################################################################
# Description: This script performs somatic variant calling using Strelka2
#              for tumor-normal paired samples. It includes configuration,
#              workflow execution, and variant normalization.
#              Strelka2 doesn't have a native tumor-only mode like Mutect2.
#
# Input:  - Recalibrated BAM files from tumor and matched normal
# Output: - Combined and normalized VCF file with somatic variants
#         - Raw SNV and indel VCF files
###############################################################################

strelka2_call() {
    # Validate required global variables
    for var in REFERENCE INTERVAL BAM STRELKA_CALL; do
        if [[ -z "${!var}" ]]; then
            echo "ERROR: Required global variable $var is not set" >&2
            return 1
        fi
    done

    local tumour=$1
    local case=$(echo "$tumour" | sed -E 's/^(.*?)-(T|N).*$/\1/')
    local normal="${case}-N"
    local out_dir="$STRELKA_CALL/${tumour}"    
    
    # Handle special case for PF samples
    if [[ "$case" =~ DFSP-PF[1-3] ]]; then
        normal=$(echo "$case" | sed 's/PF/PC/')-N
    else
        normal="${case}-N"
    fi
    
    # Define input BAM files
    local tumour_bam="${BAM}/${tumour}/${tumour}_recalibrated.bam"
    local normal_bam="${BAM}/${normal}/${normal}_recalibrated.bam"
    
    # Verify input files exist
    for bam_file in "$tumour_bam" "$normal_bam"; do
        if [[ ! -f "$bam_file" ]]; then
            echo "ERROR: BAM file not found: $bam_file" >&2
            return 1
        fi
        
        # Check if BAM index exists
        if [[ ! -f "${bam_file}.bai" ]]; then
            echo "ERROR: BAM index not found: ${bam_file}.bai" >&2
            return 1
        fi
    done
    
    # Verify reference and interval files exist
    if [[ ! -f "$REFERENCE" ]]; then
        echo "ERROR: Reference genome not found at $REFERENCE" >&2
        return 1
    fi
    
    if [[ ! -f "$INTERVAL" ]]; then
        echo "ERROR: Interval file not found at $INTERVAL" >&2
        return 1
    fi
    
    # Define indel candidates from Manta (if available)
    local indel_candidates=""
    if [[ -d "$MANTA_OUT/${case}" && -f "$MANTA_OUT/${case}/results/variants/candidateSmallIndels.vcf.gz" ]]; then
        indel_candidates="--indelCandidates $MANTA_OUT/${case}/results/variants/candidateSmallIndels.vcf.gz"
    fi
    
    # Record start time
    local start_time=$(date +%s)
    log_message "Starting Strelka2 somatic variant calling for case: $case"
    log_message "Tumor sample: $tumour"
    log_message "Normal sample: $normal"
    
    # Step 1: Configure Strelka workflow
    log_message "Configuring Strelka2 workflow..."
    
    "${STRELKA_BIN}/configureStrelkaSomaticWorkflow.py" \
        --normalBam "$normal_bam" \
        --tumorBam "$tumour_bam" \
        --referenceFasta "$REFERENCE" \
        --callRegions "$INTERVAL" \
        --exome \
        $indel_candidates \
        --runDir "$out_dir" \
        >& "${out_dir}/logs/strelka_configure.log"
    
    # Check if configuration was successful
    if [[ $? -ne 0 || ! -f "$out_dir/runWorkflow.py" ]]; then
        echo "ERROR: Strelka2 workflow configuration failed for $case" >&2
        return 1
    fi
    
    # Step 2: Run Strelka workflow
    log_message "Running Strelka2 workflow..."
    
    # Determine number of threads based on available resources
    local threads=8
    if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
        threads=$SLURM_CPUS_PER_TASK
    elif [[ -n "$PARALLEL_JOBS" ]]; then
        threads=$PARALLEL_JOBS
    fi
    
    # Run the workflow
    "$out_dir/runWorkflow.py" -m local -j $threads >& "${out_dir}/logs/strelka_workflow.log"
    
    # Check if workflow was successful
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Strelka2 workflow execution failed for $case" >&2
        return 1
    fi
    
    # Verify output files exist
    local snv_vcf="$out_dir/results/variants/somatic.snvs.vcf.gz"
    local indel_vcf="$out_dir/results/variants/somatic.indels.vcf.gz"
    
    if [[ ! -f "$snv_vcf" || ! -f "$indel_vcf" ]]; then
        echo "ERROR: Strelka2 output files not found for $case" >&2
        return 1
    fi
    
    # Step 3: Concatenate SNVs and indels
    log_message "Concatenating SNVs and indels..."
    
    local combined_vcf="$out_dir/results/variants/${case}_Strelka.vcf.gz"
    bcftools concat -Oz -a "$snv_vcf" "$indel_vcf" -o "$combined_vcf"
    
    if [[ $? -ne 0 || ! -f "$combined_vcf" ]]; then
        echo "ERROR: Failed to concatenate variant files for $case" >&2
        return 1
    fi
    
    # Step 4: Normalize variants
    log_message "Normalizing variants..."
    
    local normalized_vcf="$out_dir/results/variants/${case}_Strelka_normalized.vcf.gz"
    bcftools norm -m-both -f "$REFERENCE" -Oz -o "$normalized_vcf" "$combined_vcf"
    
    if [[ $? -ne 0 || ! -f "$normalized_vcf" ]]; then
        echo "ERROR: Failed to normalize variants for $case" >&2
        return 1
    fi
    
    # Step 5: Index the final VCF
    log_message "Indexing final VCF..."
    
    tabix "$normalized_vcf"
    
    if [[ $? -ne 0 || ! -f "${normalized_vcf}.tbi" ]]; then
        echo "ERROR: Failed to index normalized VCF for $case" >&2
        return 1
    fi
    
    # Calculate elapsed time
    local end_time=$(date +%s)
    local elapsed_time=$((end_time - start_time))
    local hours=$((elapsed_time / 3600))
    local minutes=$(( (elapsed_time % 3600) / 60 ))
    local seconds=$((elapsed_time % 60))
    
    log_message "Strelka2 variant calling completed for $case"
    log_message "Elapsed time: ${hours}h ${minutes}m ${seconds}s"
    log_message "Output: $normalized_vcf"
    
    return 0
}

export -f strelka2_call
