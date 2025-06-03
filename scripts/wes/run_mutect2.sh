#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=96:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
if ! conda activate varcall; then
    echo "ERROR: Failed to activate conda environment 'varcall'"
    exit 1
fi


# Function to run mutect2 and filtering
run_mutect2() {
    local tumour_id=$1
    local ref_dir=$2
    local bam_dir=$3
    local out_dir=$4
    local reference=$5
    local germline=$6
    local annotation_file=$7
    local interval=$8
    local pon=${9}
    
    # Extract case_id from tumour_id (everything before the second hyphen)
    case_id=$(echo $tumour_id | cut -d'-' -f1,2)
    normal_id=${case_id}-N
    
    echo "Processing sample: $tumour_id (Case: $case_id, Normal: $normal_id)"
    
    # Create output directory
    cur_dir=${out_dir}/${tumour_id}
    mkdir -p ${cur_dir}
    
    # Get Pileup Summaries
    echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries of tumour sample ..."
    gatk GetPileupSummaries \
        -I ${bam_dir}/${tumour_id}/${tumour_id}_recalibrated.bam \
        -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -O ${cur_dir}/${tumour_id}.getpileupsummaries.table

    if [[ $? -ne 0 ]]; then
        echo "ERROR: GetPileupSummaries failed for ${tumour_id}" >&2
        return 1
    fi

    #Check if presence of paired normal samples
    if [ -d ${bam_dir}/${normal_id} ]; then
        echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries of normal sample ..."
        gatk GetPileupSummaries \
            -I ${bam_dir}/${normal_id}/${normal_id}_recalibrated.bam \
            -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -O ${cur_dir}/${normal_id}.getpileupsummaries.table
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: GetPileupSummaries failed for ${normal_id}" >&2
            return 1
        fi

        # Calculate contamination on paired samples
        echo $(date +"%F") $(date +"%T") "Calculating contamination of paired samples ..."
        gatk CalculateContamination \
            -I ${cur_dir}/${tumour_id}.getpileupsummaries.table \
            -matched ${cur_dir}/${normal_id}.getpileupsummaries.table \
            -O ${cur_dir}/${tumour_id}.contamination.table \
            -segments ${cur_dir}/${tumour_id}.segments.table
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: CalculateContamination failed for ${tumour_id}_vs_${normal_id}" >&2
            return 1
        fi

        # Call variant on paired tumour samples
        echo $(date +"%F") $(date +"%T") "Calling somatic variants on paired samples ...";
        gatk --java-options -Xmx4g Mutect2 \
            -I $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            -I $bam_dir/${normal_id}/${normal_id}_recalibrated.bam \
            -normal $normal_id \
            -R ${reference} \
            -L ${interval} \
            --germline-resource ${germline}  \
            --panel-of-normals ${pon} \
            --f1r2-tar-gz ${cur_dir}/${tumour_id}.f1r2.tar.gz \
            --native-pair-hmm-threads 8 \
            --callable-depth 20 \
            -O ${cur_dir}/${tumour_id}.mutect2.vcf.gz \
            -bamout ${cur_dir}/${tumour_id}.realigned.bam \
            >& ${cur_dir}/${tumour_id}.mutect2.log

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Mutect2 failed for ${tumour_id}_vs_${normal_id}" >&2
            return 1
        fi

    else
        # Calculate contamination based on tumour samples only
        echo $(date +"%F") $(date +"%T") "Calculating contamination of tumour-only samples ..."
        gatk CalculateContamination \
            -I ${cur_dir}/${tumour_id}.getpileupsummaries.table \
            -O ${cur_dir}/${tumour_id}.contamination.table \
            -segments ${cur_dir}/${tumour_id}.segments.table
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: CalculateContamination failed for ${tumour_id}" >&2
            return 1
        fi

        # Mutect2 call variant on unpaired tumour samples
        echo $(date +"%F") $(date +"%T") "Calling somatic variants on unpaired samples ...";
        gatk --java-options -Xmx8g Mutect2 \
            -I ${bam_dir}/${tumour_id}/${tumour_id}_recalibrated.bam \
            -R ${reference} \
            -L ${interval} \
            --germline-resource ${germline}  \
            --panel-of-normals ${pon} \
            --f1r2-tar-gz ${cur_dir}/${tumour_id}.f1r2.tar.gz \
            --callable-depth 20 \
            -O ${cur_dir}/${tumour_id}.mutect2.vcf.gz \
            -bamout ${cur_dir}/${tumour_id}.realigned.bam \
            >& ${cur_dir}/${tumour_id}.mutect2.log
        
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Mutect2 failed for ${tumour_id}" >&2
            return 1
        fi
    fi

    # Learn Read Orientation Model bias for FFPE artifacts
    echo $(date +"%F") $(date +"%T") "Learning Read Orientation Model ..."
    gatk LearnReadOrientationModel \
        -I ${cur_dir}/${tumour_id}.f1r2.tar.gz \
        -O ${cur_dir}/${tumour_id}.readorientationmodel.tar.gz

    if [[ $? -ne 0 ]]; then
            echo "ERROR: LearnReadOrientationModel failed for ${tumour_id}" >&2
            return 1
    fi

    # Filter Mutect calls
    echo $(date +"%F") $(date +"%T") "Filtering Mutect calls ..."
    gatk --java-options -Xmx4g FilterMutectCalls \
        --variant ${cur_dir}/${tumour_id}.mutect2.vcf.gz \
        --stats ${cur_dir}/${tumour_id}.mutect2.vcf.gz.stats \
        --reference ${reference} \
        --ob-priors ${cur_dir}/${tumour_id}.readorientationmodel.tar.gz \
        --contamination-table ${cur_dir}/${tumour_id}.contamination.table \
        --tumor-segmentation ${cur_dir}/${tumour_id}.segments.table \
        --min-allele-fraction 0.01 \
        --unique-alt-read-count 1 \
        --output ${cur_dir}/${tumour_id}.filtermutectcalls.vcf.gz \
        >& ${cur_dir}/${tumour_id}.filtermutectcalls.log

    if [[ $? -ne 0 ]]; then
            echo "ERROR: FilterMutectCalls failed for ${tumour_id}" >&2
            return 1
    fi

    # Normalize reads to ensures variants are represented in a consistent way
    echo $(date +"%F") $(date +"%T") "Normalizing Mutect calls ..."
    bcftools norm \
        ${cur_dir}/${tumour_id}.filtermutectcalls.vcf.gz \
        -m-both -f ${reference} \
        -Oz -o ${cur_dir}/${tumour_id}.normalized.vcf.gz 

    # keep variants that have passed all filters (low-quality or failed variant calls)
    echo $(date +"%F") $(date +"%T") "Filtering PASS variants ..."
    bcftools view \
        -f PASS ${cur_dir}/${tumour_id}.normalized.vcf.gz \
        -o ${cur_dir}/${tumour_id}.passed.vcf.gz

    # Annotate repeatmasker and blacklist regions
    echo $(date +"%F") $(date +"%T") "Annotating repeatmasker regions ..."
    bcftools annotate \
        ${cur_dir}/${tumour_id}.passed.vcf.gz \
        --header-lines ${out_dir}/vcf.rm.header \
        --annotations ${ref_dir}/RepeatMasker.bed.gz \
        --columns CHROM,FROM,TO,RepeatMasker \
        --output ${cur_dir}/${tumour_id}.repeatmasker.vcf.gz

    echo $(date +"%F") $(date +"%T") "Annotating blacklist regions ..."
    bcftools annotate \
        ${cur_dir}/${tumour_id}.repeatmasker.vcf.gz \
        --header-lines ${out_dir}/vcf.map.header \
        --annotations ${ref_dir}/blacklist.bed.gz \
        --columns CHROM,FROM,TO,EncodeDacMapability \
        --output-type z \
        --output ${cur_dir}/${tumour_id}.blacklist.vcf.gz
    
    # Filter out variants in RepeatMasker or Mapability
    echo $(date +"%F") $(date +"%T") "Filtering RepeatMasker and blacklist regions ..."
    bcftools filter \
        ${cur_dir}/${tumour_id}.blacklist.vcf.gz \
        -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \
        -Oz \
        -o ${cur_dir}/${tumour_id}.final.vcf.gz

    # Index the final VCF file
    tabix -p vcf ${cur_dir}/${tumour_id}.final.vcf.gz

    # Annotate variants by Funcotator
    echo $(date +"%F") $(date +"%T") "Annotating variants with Funcotator ..."
    gatk Funcotator \
        -R ${reference} \
        -V ${cur_dir}/${tumour_id}.final.vcf.gz \
        -O ${cur_dir}/${tumour_id}.funcotator.maf.gz \
        -L ${interval} \
        --output-file-format MAF \
        --data-sources-path ${annotation_file} \
        --ref-version hg38 \
        --remove-filtered-variants true \
        >& ${cur_dir}/${tumour_id}.funcotator.log

    if [[ $? -ne 0 ]]; then
            echo "ERROR: Funcotator failed for ${tumour_id}" >&2
            return 1
    fi

    # Extract relevant columns from the Funcotator output
    echo $(date +"%F") $(date +"%T") "Extracting Funcotator data ..."
    less -S ${cur_dir}/${tumour_id}.funcotator.maf.gz | \
        grep -v "#" > ${cur_dir}/${tumour_id}.funcotator.tsv

    # Annotate variant by Annovar
    echo $(date +"%F") $(date +"%T") "Annotating variants with Annovar ..."
    perl ${ref_dir}/annovar/table_annovar.pl \
        ${cur_dir}/${tumour_id}.final.vcf.gz \
        ${ref_dir}/annovar/humandb/ \
        -buildver hg38 \
        -out ${cur_dir}/${tumour_id} \
        -remove \
        -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
        -operation gx,r,f,f,f,f,f \
        -nastring . \
        -polish \
        -xreffile ${ref_dir}/annovar/example/gene_fullxref.txt \
        --otherinfo \
        --vcfinput \
        >& ${cur_dir}/${tumour_id}.annovar.log

    if [[ $? -ne 0 ]]; then
            echo "ERROR: Annovar failed for ${tumour_id}" >&2
            return 1
    fi

    # Extract relevant columns from the Annovar output
    echo $(date +"%F") $(date +"%T") "Processing Annovar output ..."
    less -S ${cur_dir}/${tumour_id}.hg38_multianno.txt | \
        awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > \
        ${cur_dir}/${tumour_id}.annovar.txt
    
    # Clean up intermediate files
    echo $(date +"%F") $(date +"%T") "Cleaning up intermediate files ..."
    rm -f $cur_dir/${tumour_id}.normalized.vcf.gz
    rm -f $cur_dir/${tumour_id}.passed.vcf.gz
    rm -f $cur_dir/${tumour_id}.repeatmasker.vcf.gz
    rm -f $cur_dir/${tumour_id}.blacklist.vcf.gz

    echo $(date +"%F") $(date +"%T") "Completed processing sample: $tumour_id"
}

# Export function to make it available to GNU parallel
export -f run_mutect2

# Define directories
export ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
export bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
export out_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_tumour_normal"
mkdir -p ${out_dir}

# Define reference files
export reference=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
export germline=${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
export annotation_file=${ref_dir}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
export interval=${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
export pon=${ref_dir}/pon_dfsp/pon.vcf.gz

echo "Reference directory:  ${ref_dir}"
echo "Bam directory:        ${bam_dir}"
echo "Work directory:       ${out_dir}"
echo "Reference:            ${reference}"
echo "Germline:             ${germline}"
echo "Annotation file:      ${annotation_file}"
echo "Interval:             ${interval}"
echo "Panel of Normal:      ${pon}"

# Create header files needed for repeatmasker and blacklist annotation
echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > ${out_dir}/vcf.rm.header
echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > ${out_dir}/vcf.map.header

# Get list of tumor samples to process
#vcf_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2"
# find "$vcf_dir" -name "*_unfiltered.vcf.gz" | sort | sed 's|.*/||' | sed 's/_unfiltered.vcf.gz$//' > "${out_dir}/sample_list.txt"
tumour_samples="/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/tumour_all_samples.txt"

# Count the number of samples
num_sample=$(cat "$tumour_samples" | wc -l)
echo "Number of samples = $num_sample"

# Calculate number of parallel jobs
if [ $num_sample -ge 15 ]; then
    jobs=15
else
    jobs=$num_sample
fi

# Run the processing in parallel
cat "${out_dir}/sample_list.txt" | parallel \
    --jobs ${jobs} \
    --progress \
    run_mutect2 {} "$ref_dir" "$bam_dir" "$out_dir" "$reference" "$germline" "$annotation_file" "$interval" "$pon"

rm ${out_dir}/vcf.rm.header
rm ${out_dir}/vcf.map.header

echo "All samples processed successfully."