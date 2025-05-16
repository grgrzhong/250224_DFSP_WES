#!/bin/bash

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

# Define directories
export ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
export bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
export vcf_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2"
export work_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_tumour_normal"
mkdir -p ${work_dir}

# Define reference files
export REFERENCE=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
export GERMLINE=${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
export ANNOTATION_FILE=${ref_dir}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
export INTERVAL=${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
export PON=${ref_dir}/pon_dfsp/pon.vcf.gz

echo "Reference directory:  ${ref_dir}"
echo "Bam directory:        ${bam_dir}"
echo "Work directory:       ${work_dir}"
echo "Reference:            ${REFERENCE}"
echo "Germline:             ${GERMLINE}"
echo "Annotation file:      ${ANNOTATION_FILE}"
echo "Interval:             ${INTERVAL}"
echo "PON:                  ${PON}"


# Function to run mutect2 and filtering
mutect_call_filter() {
    local tumour_id=$1
    local ref_dir=$2
    local bam_dir=$3
    local vcf_dir=$4
    local work_dir=$5
    local REFERENCE=$6
    local GERMLINE=$7
    local ANNOTATION_FILE=$8
    local INTERVAL=$9
    local PON=${10}
    
    # Extract case_id from tumour_id (everything before the second hyphen)
    case_id=$(echo $tumour_id | cut -d'-' -f1,2)
    normal_id=${case_id}-N
    
    echo "Processing sample: $tumour_id (Case: $case_id, Normal: $normal_id)"
    
    # Create output directory
    mkdir -p ${work_dir}/${tumour_id}
    VAROUT=${work_dir}/${tumour_id}
    
    # Get Pileup Summaries
    echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries of tumour sample ..."
    gatk GetPileupSummaries \
        -I $bam_dir/$tumour_id/${tumour_id}_recalibrated.bam \
        -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -O ${VAROUT}/${tumour_id}.getpileupsummaries.table

    #Check if presence of paired normal samples
    if [ -d ${bam_dir}/${normal_id} ]; then
        echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries of normal sample ..."
        gatk GetPileupSummaries \
            -I $bam_dir/$normal_id/${normal_id}_recalibrated.bam \
            -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -O ${VAROUT}/${normal_id}.getpileupsummaries.table

        # Calculate contamination on paired samples
        echo $(date +"%F") $(date +"%T") "Calculating contamination of paired samples ..."
        gatk CalculateContamination \
            -I ${VAROUT}/${tumour_id}.getpileupsummaries.table \
            -matched ${VAROUT}/${normal_id}.getpileupsummaries.table \
            -O ${VAROUT}/${tumour_id}.contamination.table \
            -segments ${VAROUT}/${tumour_id}.segments.table

        # Call variant on paired tumour samples
        echo $(date +"%F") $(date +"%T") "Calling somatic variants on paired samples ...";
        gatk --java-options -Xmx8g Mutect2 \
            -I $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            -I $bam_dir/${normal_id}/${normal_id}_recalibrated.bam \
            -normal $normal_id \
            -R $REFERENCE \
            -L $INTERVAL \
            --germline-resource $GERMLINE  \
            --panel-of-normals $PON \
            --f1r2-tar-gz ${VAROUT}/${tumour_id}.f1r2.tar.gz \
            --native-pair-hmm-threads 8 \
            --callable-depth 20 \
            -O ${VAROUT}/${tumour_id}.mutect2.vcf.gz \
            -bamout ${VAROUT}/${tumour_id}.realigned.bam \
            >& ${VAROUT}/${tumour_id}.mutect2.log
    else
        # Calculate contamination based on tumour samples only
        echo $(date +"%F") $(date +"%T") "Calculating contamination of tumour-only samples ..."
        gatk CalculateContamination \
            -I ${VAROUT}/${tumour_id}.getpileupsummaries.table \
            -O ${VAROUT}/${tumour_id}.contamination.table \
            -segments ${VAROUT}/${tumour_id}.segments.table
        
        # Mutect2 call variant on unpaired tumour samples
        echo $(date +"%F") $(date +"%T") "Calling somatic variants on unpaired samples ...";
        gatk --java-options -Xmx8g Mutect2 \
            -I $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            -R $REFERENCE \
            -L $INTERVAL \
            --germline-resource $GERMLINE  \
            --panel-of-normals $PON \
            --f1r2-tar-gz ${VAROUT}/${tumour_id}.f1r2.tar.gz \
            --callable-depth 20 \
            -O ${VAROUT}/${tumour_id}.mutect2.vcf.gz \
            -bamout ${VAROUT}/${tumour_id}.realigned.bam \
            >& ${VAROUT}/${tumour_id}.mutect2.log
    fi

    # Learn Read Orientation Model bias for FFPE artifacts
    gatk LearnReadOrientationModel \
        -I ${VAROUT}/${tumour_id}/${tumour_id}.f1r2.tar.gz \
        -O ${VAROUT}/${tumour_id}.readorientationmodel.tar.gz

    # Filter Mutect calls
    echo $(date +"%F") $(date +"%T") "Filtering Mutect calls ..."
    gatk --java-options -Xmx4g FilterMutectCalls \
        --variant ${VAROUT}/${tumour_id}/${tumour_id}.mutect2.vcf.gz \
        --stats ${VAROUT}/${tumour_id}/${tumour_id}.mutect2.vcf.gz.stats \
        --reference $REFERENCE \
        --ob-priors $VAROUT/${tumour_id}.readorientationmodel.tar.gz \
        --contamination-table $VAROUT/${tumour_id}.contamination.table \
        --tumor-segmentation $VAROUT/${tumour_id}.segments.table \
        --min-allele-fraction 0.01 \
        --unique-alt-read-count 1 \
        --output $VAROUT/${tumour_id}.filtermutectcalls.vcf.gz \
        >& ${VAROUT}/${tumour_id}.filtermutectcalls.log

    # Normalize reads to ensures variants are represented in a consistent way
    echo $(date +"%F") $(date +"%T") "Normalizing Mutect calls ..."
    bcftools norm \
        ${VAROUT}/${tumour_id}.filtermutectcalls.vcf.gz \
        -m-both -f $REFERENCE \
        -Oz -o ${VAROUT}/${tumour_id}.normalized.vcf.gz 

    # keep variants that have passed all filters (low-quality or failed variant calls)
    echo $(date +"%F") $(date +"%T") "Normalizing filtered Mutect calls ..."
    bcftools view \
        -f PASS ${VAROUT}/${tumour_id}.normalized.vcf.gz \
        -o ${VAROUT}/${tumour_id}.passed.vcf.gz

    # Annotate repeatmasker and blacklist regions
    echo $(date +"%F") $(date +"%T") "Annotating repeatmasker regions ..."
    bcftools annotate \
        ${VAROUT}/${tumour_id}.passed.vcf.gz \
        --header-lines ${work_dir}/vcf.rm.header \
        --annotations ${ref_dir}/RepeatMasker.bed.gz \
        --columns CHROM,FROM,TO,RepeatMasker \
        --output ${VAROUT}/${tumour_id}.repeatmasker.vcf.gz

    echo $(date +"%F") $(date +"%T") "Annotating blacklist regions ..."
    bcftools annotate \
        ${VAROUT}/${tumour_id}.repeatmasker.vcf.gz \
        --header-lines ${work_dir}/vcf.map.header \
        --annotations ${ref_dir}/blacklist.bed.gz \
        --columns CHROM,FROM,TO,EncodeDacMapability \
        --output-type z \
        --output ${VAROUT}/${tumour_id}.blacklist.vcf.gz
    
    # Filter out variants in RepeatMasker or Mapability
    echo $(date +"%F") $(date +"%T") "Filtering RepeatMasker and blacklist regions ..."
    bcftools filter \
        ${VAROUT}/${tumour_id}.blacklist.vcf.gz \
        -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \
        -Oz \
        -o ${VAROUT}/${tumour_id}.final.vcf.gz

    tabix ${VAROUT}/${tumour_id}.final.vcf.gz

    # Annotate variants by Funcotator
    echo $(date +"%F") $(date +"%T") "Annotating variants with Funcotator ..."
    gatk Funcotator \
        -R $REFERENCE \
        -V ${VAROUT}/${tumour_id}.final.vcf.gz \
        -O ${VAROUT}/${tumour_id}.funcotator.maf.gz \
        -L $INTERVAL \
        --output-file-format MAF \
        --data-sources-path $ANNOTATION_FILE \
        --ref-version hg38 \
        --remove-filtered-variants true \
        >& ${VAROUT}/${tumour_id}.funcotator.log
    
    # Extract relevant columns from the Funcotator output
    less -S ${VAROUT}/${tumour_id}.funcotator.maf.gz | \
        grep -v "#" > ${VAROUT}/${tumour_id}.funcotator.tsv

    # Annotate variant by Annovar
    echo $(date +"%F") $(date +"%T") "Annotating variants with Annovar ..."
    perl ${ref_dir}/annovar/table_annovar.pl \
        ${VAROUT}/${tumour_id}.final.vcf.gz \
        ${ref_dir}/annovar/humandb/ \
        -buildver hg38 \
        -out ${VAROUT}/${tumour_id} \
        -remove \
        -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
        -operation gx,r,f,f,f,f,f \
        -nastring . \
        -polish \
        -xreffile ${ref_dir}/annovar/example/gene_fullxref.txt \
        --otherinfo \
        --vcfinput \
        >& ${VAROUT}/${tumour_id}.annovar.log

    # Extract relevant columns from the Annovar output
    less -S ${VAROUT}/${tumour_id}.hg38_multianno.txt | \
        awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > \
        ${VAROUT}/${tumour_id}.annovar.txt
    
    # Clean up intermediate files
    # rm $VAROUT/${tumour_id}.filtermutectcalls.vcf.gz
    # rm $VAROUT/${tumour_id}.filtermutectcalls.vcf.gz.tbi
    rm $VAROUT/${tumour_id}.normalized.vcf.gz
    rm $VAROUT/${tumour_id}.passed.vcf.gz
    rm $VAROUT/${tumour_id}.repeatmasker.vcf.gz
    rm $VAROUT/${tumour_id}.blacklist.vcf.gz

    echo $(date +"%F") $(date +"%T") "Completed processing sample: $tumour_id"
}

# Export function to make it available to GNU parallel
export -f mutect_call_filter

# Create header files needed for repeatmasker and blacklist annotation
echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > ${work_dir}/vcf.rm.header
echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > ${work_dir}/vcf.map.header

# Get list of tumor samples to process
# You can customize this part based on how you want to identify samples
# Example 1: Specific list
# tumour_ids="DFSP-028-T DFSP-029-T DFSP-030-T-P1"

# Example 2: Find all tumor samples from the vcf directory
find "$vcf_dir" -name "*_unfiltered.vcf.gz" | sort | sed 's|.*/||' | sed 's/_unfiltered.vcf.gz$//' > "${work_dir}/sample_list.txt"
# find "$bam_dir" -name "*_recalibrated.bam" | sort | sed 's|.*/||' | sed 's/_recalibrated.bam$//' > "${work_dir}/sample_list.txt"
# cat "${work_dir}/sample_list.txt"

# Number of parallel processes to run (adjust based on your system's capacity)
PARALLEL_JOBS=30

# Run the processing in parallel
echo "Starting parallel processing of samples with $PARALLEL_JOBS jobs..."
cat "${work_dir}/sample_list.txt" | parallel \
    --jobs $PARALLEL_JOBS \
    --progress \
    --joblog ${work_dir}/parallel.log \
    mutect_call_filter {} "$ref_dir" "$bam_dir" "$vcf_dir" "$work_dir" "$REFERENCE" "$GERMLINE" "$ANNOTATION_FILE" "$INTERVAL" "$PON"

echo "All samples processed successfully."
rm ${work_dir}/sample_list.txt
rm ${work_dir}/vcf.rm.header
rm ${work_dir}/vcf.map.header