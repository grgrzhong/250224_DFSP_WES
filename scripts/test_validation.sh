#!/bin/bash

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

# Define directories
export ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
export bam_dir="/path/to/NA12878/bam" # Update with NA12878 BAM location
export vcf_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2"
export work_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/validation/NA12878"
export truth_vcf="/path/to/NA12878/giab/HG001_GRCh38_GIAB_highconf.vcf.gz" # GIAB truth set
export truth_bed="/path/to/NA12878/giab/HG001_GRCh38_GIAB_highconf.bed" # GIAB confidence regions
mkdir -p ${work_dir}

# Define reference files (same as your existing pipeline)
export REFERENCE=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
export GERMLINE=${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
export ANNOTATION_FILE=${ref_dir}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
export INTERVAL=${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
export PON=${ref_dir}/pon_dfsp/pon.vcf.gz

# Set NA12878 sample ID 
export na12878_id="NA12878"

# Function for NA12878 validation
validate_na12878() {
    local na12878_id="NA12878"
    local ref_dir=$1
    local bam_dir=$2
    local work_dir=$3
    local REFERENCE=$4
    local GERMLINE=$5
    local ANNOTATION_FILE=$6
    local INTERVAL=$7
    local truth_vcf=$8
    local truth_bed=$9
    
    echo "Validating Mutect2 pipeline with NA12878 reference sample"
    
    # Create output directory
    mkdir -p ${work_dir}/${na12878_id}
    VALOUT=${work_dir}/${na12878_id}
    
    # Step 1: GetPileupSummaries
    echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries for NA12878..."
    gatk GetPileupSummaries \
        -I $bam_dir/${na12878_id}/${na12878_id}_recalibrated.bam \
        -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -O ${VALOUT}/${na12878_id}.getpileupsummaries.table

    # Step 2: Calculate Contamination (single sample mode)
    echo $(date +"%F") $(date +"%T") "Calculating contamination for NA12878..."
    gatk CalculateContamination \
        -I ${VALOUT}/${na12878_id}.getpileupsummaries.table \
        -O ${VALOUT}/${na12878_id}.contamination.table \
        -segments ${VALOUT}/${na12878_id}.segments.table
    
    # Step 3: Run Mutect2 (using NA12878 as "tumor" sample)
    echo $(date +"%F") $(date +"%T") "Running Mutect2 on NA12878..."
    gatk --java-options -Xmx8g Mutect2 \
        -I $bam_dir/${na12878_id}/${na12878_id}_recalibrated.bam \
        -R $REFERENCE \
        -L $INTERVAL \
        --germline-resource $GERMLINE \
        --f1r2-tar-gz ${VALOUT}/${na12878_id}.f1r2.tar.gz \
        --callable-depth 20 \
        -O ${VALOUT}/${na12878_id}_unfiltered.vcf.gz \
        >& ${VALOUT}/${na12878_id}.mutect2.log
    
    # Step 4: Learn Read Orientation Model
    echo $(date +"%F") $(date +"%T") "Learning read orientation model..."
    gatk LearnReadOrientationModel \
        -I ${VALOUT}/${na12878_id}.f1r2.tar.gz \
        -O ${VALOUT}/${na12878_id}.read-orientation-model.tar.gz
    
    # Step 5: Filter Mutect Calls
    echo $(date +"%F") $(date +"%T") "Filtering Mutect calls..."
    gatk --java-options -Xmx4g FilterMutectCalls \
        --variant ${VALOUT}/${na12878_id}_unfiltered.vcf.gz \
        --stats ${VALOUT}/${na12878_id}_unfiltered.vcf.gz.stats \
        --reference $REFERENCE \
        --ob-priors ${VALOUT}/${na12878_id}.read-orientation-model.tar.gz \
        --contamination-table ${VALOUT}/${na12878_id}.contamination.table \
        --tumor-segmentation ${VALOUT}/${na12878_id}.segments.table \
        --min-allele-fraction 0.01 \
        --unique-alt-read-count 1 \
        --output ${VALOUT}/${na12878_id}_filtered.vcf.gz \
        >& ${VALOUT}/${na12878_id}.filtermutectcalls.log
    
    # Continue with the rest of your pipeline (normalization, annotation, etc.)
    # ...
    
    # Add validation step - compare against truth set
    echo $(date +"%F") $(date +"%T") "Comparing results against truth set..."
    hap.py \
        $truth_vcf \
        ${VALOUT}/${na12878_id}.final.vcf.gz \
        -f $truth_bed \
        -r $REFERENCE \
        -o ${VALOUT}/${na12878_id}_happy \
        --engine=vcfeval \
        --pass-only \
        --target-regions=$INTERVAL
    
    # Output validation summary
    echo "Validation Results Summary:"
    echo "=========================="
    cat ${VALOUT}/${na12878_id}_happy.summary.csv
}