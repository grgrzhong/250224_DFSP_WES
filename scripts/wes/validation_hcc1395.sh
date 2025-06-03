#!/bin/bash

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

# Define reference files
export ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
export dbsnp="${ref_dir}/Population_database/dbSNP.vcf.gz"
export pon=${ref_dir}/pon_dfsp/pon.vcf.gz
# export interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
# export bait_interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list"
# export target_interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list"
export germline="${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz"

# export interval="/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395/raw/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed.gz"
export reference="/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395/raw/GRCh38/GRCh38.d1.vd1.fa"

sample=HCC1395
tumour_id=WES_IL_T_1
normal_id=WES_IL_N_1

work_dir=/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/${sample}
bam_dir="${work_dir}/raw"

# fastq_1=/home/zhonggr/projects/250224_DFSP_WES/data/reference/benchmark/NA12878/raw/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
# fastq_2=/home/zhonggr/projects/250224_DFSP_WES/data/reference/benchmark/NA12878/raw/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# zcat ${fastq_1} | head -n 15
# zcat ${fastq_2} | head -n 15
# fastq_trimmed_dir="${work_dir}/preprocessing/fastq_trimmed"
# mkdir -p ${fastq_trimmed_dir}
# bam_dir="${work_dir}/preprocessing/bam"
# mkdir -p ${bam_dir}
depth=20
mutect2_dir="${work_dir}/mutect2_${depth}"
mkdir -p ${mutect2_dir}

echo "reference: ${reference}"
echo "dbsnp:     ${dbsnp}"
echo "pon:       ${pon}"
echo "germline:  ${germline}"
echo "tumour_id: ${tumour_id}"
echo "normal_id: ${normal_id}"
echo "work_dir:  ${work_dir}"
echo "mutect2_dir: ${mutect2_dir}"

# # Experimental UMI options
# umi_opts=""
# use_umi=false
# umi_loc="per_read"
# umi_len=8
# if [[ "${use_umi}" == "true" ]]; then
#     umi_opts="--umi --umi_loc ${umi_loc} --umi_len ${umi_len}"
# fi


# # adapter_1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
# # adapter_2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# adapter_1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
# adapter_2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"

# ############################################################################
# ## Preprocessing and Alignment ----
# ############################################################################
# # Step 1: Fastp trimming
# echo $(date +"%F") $(date +"%T") "Running Fastp trimming for $sample..."
# fastp \
#     -i ${fastq_1} \
#     -I ${fastq_2} \
#     -o ${fastq_trimmed_dir}/${sample}_trimmed_1.fastq.gz \
#     -O ${fastq_trimmed_dir}/${sample}_trimmed_2.fastq.gz \
#     --detect_adapter_for_pe \
#     --adapter_sequence=${adapter_1} \
#     --adapter_sequence_r2=${adapter_2} \
#     --trim_poly_g \
#     --trim_poly_x \
#     ${umi_opts} \
#     -j ${fastq_trimmed_dir}/${sample}.json \
#     -h ${fastq_trimmed_dir}/${sample}.html \
#     -w 8

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: Fastp trimming failed for ${sample}" >&2
#     return 1
# fi

# # Step 2: BWA Alignment
# echo $(date +"%F") $(date +"%T") "Aligning ${sample} to reference genome..."

# kit_name="Nextera_Rapid_Capture_Exome"
# platform="ILLUMINA"
# platform_model="HISEQ2500" 

# bwa mem -M -t 16 \
#     -R "@RG\tID:${sample}\tLB:${kit_name}\tPL:${platform}\tPM:${platform_model}\tSM:${sample}\tPU:NA" \
#     ${reference} \
#     ${fastq_trimmed_dir}/${sample}_trimmed_1.fastq.gz \
#     ${fastq_trimmed_dir}/${sample}_trimmed_2.fastq.gz > ${bam_dir}/${sample}.sam

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: BWA alignment failed for ${sample}" >&2
#     return 1
# fi

# samtools view -Sb ${bam_dir}/${sample}.sam > ${bam_dir}/${sample}.bwa.bam

# # # Step 3: Extract and tag UMI, remove as 
# # echo $(date +"%F") $(date +"%T") "Extract and tag UMI for ${sample}..."
# # python /home/zhonggr/projects/250224_DFSP_WES/bin/python/tag_umi.py \
# #     -I ${bam_dir}/${sample}.bwa.bam \
# #     -O ${bam_dir}/${sample}.umi.bam

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: UMI tagging failed for ${sample}" >&2
#     return 1
# fi

# # Step 4: Sort by coordinate
# echo $(date +"%F") $(date +"%T") "Sorting BAM file for ${sample}..."
# samtools sort \
#     ${bam_dir}/${sample}.bwa.bam \
#     -o ${bam_dir}/${sample}.sorted.bam

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: BAM sorting failed for ${sample}" >&2
#     return 1
# fi

# # Step 5: Mark duplicates
# echo $(date +"%F") $(date +"%T") "Marking duplicates for ${sample}..."
# gatk --java-options -Xmx4g MarkDuplicates \
#     -I ${bam_dir}/${sample}.sorted.bam \
#     -M ${bam_dir}/${sample}.metrics.txt \
#     -O ${bam_dir}/${sample}.marked.bam \
#     >& ${bam_dir}/${sample}.markduplicates.log
#     # Remove: --BARCODE_TAG "RX"  # No UMIs available

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: Marking duplicates failed for ${sample}" >&2
#     return 1
# fi

# # Step 6: Index BAM
# echo $(date +"%F") $(date +"%T") "Indexing BAM file for ${sample}..."
# samtools index ${bam_dir}/${sample}.marked.bam

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: BAM indexing failed for ${sample}" >&2
#     return 1
# fi

# # Step 7: Base recalibration
# echo $(date +"%F") $(date +"%T") "Running base recalibration for ${sample}..."
# gatk BaseRecalibrator \
#     -I ${bam_dir}/${sample}.marked.bam \
#     -R ${reference} \
#     -L ${interval} \
#     -O ${bam_dir}/${sample}.recal.table \
#     --known-sites ${dbsnp} \
#     >& ${bam_dir}/${sample}.baserecalibrator.log

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: Base recalibration failed for ${sample}" >&2
#     return 1
# fi

# # Step 8: Apply BQSR
# echo $(date +"%F") $(date +"%T") "Applying BQSR for ${sample}..."
# gatk ApplyBQSR \
#     -I ${bam_dir}/${sample}.marked.bam \
#     -O ${bam_dir}/${sample}.recal.bam \
#     -L ${interval} \
#     -bqsr ${bam_dir}/${sample}.recal.table \
#     --create-output-bam-md5 \
#     >& ${bam_dir}/${sample}.applybqsr.log

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: Applying BQSR failed for ${sample}" >&2
#     return 1
# fi

# # Step 9: Collect metrics
# echo $(date +"%F") $(date +"%T") "Collecting HsMetrics for ${sample}..."
# gatk CollectHsMetrics \
#     -I ${bam_dir}/${sample}.recal.bam \
#     -O ${bam_dir}/${sample}.hsmetrics.txt \
#     -R ${reference} \
#     -BI ${bait_interval} \
#     -TI ${target_interval} \
#     >& ${bam_dir}/${sample}.collecthsmetrics.log

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: Collecting HsMetrics failed for ${sample}" >&2
#     return 1
# fi

# # Step 10: Generate alignment stats
# echo $(date +"%F") $(date +"%T") "Generating alignment stats for ${sample}..."
# bamtools stats \
#     -in ${bam_dir}/${sample}.recal.bam \
#     > ${bam_dir}/${sample}.alnstat.txt

# if [[ $? -ne 0 ]]; then
#     echo "ERROR: Generating alignment stats failed for ${sample}" >&2
#     return 1
# fi

# # Clean up intermediate files
# echo $(date +"%F") $(date +"%T") "Cleaning up intermediate files for ${sample}..."
# rm -f \
#     ${bam_dir}/${sample}.bwa.bam \
#     ${bam_dir}/${sample}.umi.bam \
#     ${bam_dir}/${sample}.sorted.bam \
#     ${bam_dir}/${sample}.marked.bam \
#     ${bam_dir}/${sample}.marked.bam.bai

# echo $(date +"%F") $(date +"%T") "Completed processing sample: ${sample}"

#############################################################################
## Mutect2 calling
#############################################################################
# Create header files needed for repeatmasker and blacklist annotation
echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > ${mutect2_dir}/vcf.rm.header
echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > ${mutect2_dir}/vcf.map.header

# # Step 1: GetPileupSummaries
echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries ..."
gatk GetPileupSummaries \
    -I ${work_dir}/raw/${tumour_id}.bwa.dedup.bam \
    -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
    -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
    -O ${mutect2_dir}/${tumour_id}.getpileupsummaries.table \
    >& ${mutect2_dir}/${tumour_id}.getpileupsummaries.log

gatk GetPileupSummaries \
    -I ${work_dir}/raw/${normal_id}.bwa.dedup.bam \
    -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
    -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
    -O ${mutect2_dir}/${normal_id}.getpileupsummaries.table \
    >& ${mutect2_dir}/${normal_id}.getpileupsummaries.log

# Step 2: Calculate Contamination
echo $(date +"%F") $(date +"%T") "Calculating contamination for NA12878..."
gatk CalculateContamination \
    -I ${mutect2_dir}/${tumour_id}.getpileupsummaries.table \
    -matched ${mutect2_dir}/${normal_id}.getpileupsummaries.table \
    -O ${mutect2_dir}/${tumour_id}.contamination.table \
    -segments ${mutect2_dir}/${tumour_id}.segments.table \
    >& ${mutect2_dir}/${tumour_id}.calculatecontamination.log
    
# Step 3: Run Mutect2
echo $(date +"%F") $(date +"%T") "Running Mutect2 calling for ${sample}..."
gatk --java-options -Xmx8g Mutect2 \
    -I ${work_dir}/raw/${tumour_id}.bwa.dedup.bam \
    -I ${work_dir}/raw/${normal_id}.bwa.dedup.bam \
    -normal ${normal_id} \
    -R ${reference} \
    --germline-resource ${germline} \
    --panel-of-normals ${pon} \
    --f1r2-tar-gz ${mutect2_dir}/${tumour_id}.f1r2.tar.gz \
    --native-pair-hmm-threads 8 \
    --callable-depth ${depth} \
    -O ${mutect2_dir}/${tumour_id}.mutect2.vcf.gz \
    >& ${mutect2_dir}/${tumour_id}.mutect2.log

raw_variants=$(bcftools view -H ${mutect2_dir}/${tumour_id}.mutect2.vcf.gz | wc -l)
echo "Raw variants called by Mutect2: $raw_variants"

# Step 4: Learn Read Orientation Model
echo $(date +"%F") $(date +"%T") "Learning read orientation model..."
gatk LearnReadOrientationModel \
    -I ${mutect2_dir}/${tumour_id}.f1r2.tar.gz \
    -O ${mutect2_dir}/${tumour_id}.readorientationmodel.tar.gz \
    >& ${mutect2_dir}/${tumour_id}.learnreadorientationmodel.log
    
# Step 5: Filter Mutect Calls
echo $(date +"%F") $(date +"%T") "Filtering Mutect calls..."
gatk --java-options -Xmx4g FilterMutectCalls \
    --variant ${mutect2_dir}/${tumour_id}.mutect2.vcf.gz \
    --stats ${mutect2_dir}/${tumour_id}.mutect2.vcf.gz.stats \
    --reference ${reference} \
    --ob-priors ${mutect2_dir}/${tumour_id}.readorientationmodel.tar.gz \
    --contamination-table ${mutect2_dir}/${tumour_id}.contamination.table \
    --tumor-segmentation ${mutect2_dir}/${tumour_id}.segments.table \
    --min-allele-fraction 0.01 \
    --unique-alt-read-count 1 \
    --output ${mutect2_dir}/${tumour_id}.filtermutectcalls.vcf.gz \
    >& ${mutect2_dir}/${tumour_id}.filtermutectcalls.log

filtered_variants=$(bcftools view -H ${mutect2_dir}/${tumour_id}.filtermutectcalls.vcf.gz | wc -l)
echo "After FilterMutectCalls: $filtered_variants (removed: $((raw_variants - filtered_variants)))"

# # Normalize reads to ensures variants are represented in a consistent way
echo $(date +"%F") $(date +"%T") "Normalizing Mutect calls ..."
bcftools norm \
    ${mutect2_dir}/${tumour_id}.filtermutectcalls.vcf.gz \
    -m-both -f ${reference} \
    -Oz -o ${mutect2_dir}/${tumour_id}.normalized.vcf.gz 

tabix -p vcf ${mutect2_dir}/${tumour_id}.normalized.vcf.gz

norm_variants=$(bcftools view -H ${mutect2_dir}/${tumour_id}.normalized.vcf.gz | wc -l)
echo "After normalize: $norm_variants (removed: $((filtered_variants - norm_variants)))"

# keep variants that have passed all filters (low-quality or failed variant calls)
echo $(date +"%F") $(date +"%T") "Filtering PASS variants ..."
bcftools view \
    -f PASS ${mutect2_dir}/${tumour_id}.normalized.vcf.gz \
    -o ${mutect2_dir}/${tumour_id}.passed.vcf.gz

pass_variants=$(bcftools view -H ${mutect2_dir}/${tumour_id}.passed.vcf.gz | wc -l)
echo "After pass filtering: $pass_variants (removed: $((norm_variants - pass_variants)))"

# # Annotate repeatmasker and blacklist regions
# echo $(date +"%F") $(date +"%T") "Annotating repeatmasker regions ..."
bcftools annotate \
    ${mutect2_dir}/${tumour_id}.passed.vcf.gz \
    --header-lines ${mutect2_dir}/vcf.rm.header \
    --annotations ${ref_dir}/RepeatMasker.bed.gz \
    --columns CHROM,FROM,TO,RepeatMasker \
    --output ${mutect2_dir}/${tumour_id}.repeatmasker.vcf.gz

echo $(date +"%F") $(date +"%T") "Annotating blacklist regions ..."
bcftools annotate \
    ${mutect2_dir}/${tumour_id}.repeatmasker.vcf.gz \
    --header-lines ${mutect2_dir}/vcf.map.header \
    --annotations ${ref_dir}/blacklist.bed.gz \
    --columns CHROM,FROM,TO,EncodeDacMapability \
    --output-type z \
    --output ${mutect2_dir}/${tumour_id}.blacklist.vcf.gz

# # Filter out variants in RepeatMasker or Mapability
echo $(date +"%F") $(date +"%T") "Filtering RepeatMasker and blacklist regions ..."
bcftools filter \
    ${mutect2_dir}/${tumour_id}.blacklist.vcf.gz \
    -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \
    -Oz \
    -o ${mutect2_dir}/${tumour_id}.final.vcf.gz

# Index the final VCF file
tabix -p vcf ${mutect2_dir}/${tumour_id}.final.vcf.gz

final_variants=$(bcftools view -H ${mutect2_dir}/${tumour_id}.final.vcf.gz | wc -l)
echo "After pass filtering: $final_variants (removed: $((pass_variants - final_variants)))"

###############################################################################
## Hap.py the mutect2 call vs  truth set
###############################################################################
echo $(date +"%F") $(date +"%T") "Comparing results against truth set..."

## The original mutect call variants=20235
truth_mutect2_vcf=/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395/raw/WES_IL_1.novo.muTect2.vcf.gz
truth_mutect2_variants=$(bcftools view -H ${truth_mutect2_vcf} | wc -l)
echo "Mutect2 variants: $truth_mutect2_variants"
bcftools query -l ${truth_mutect2_vcf}
zgrep "^#CHROM" ${truth_mutect2_vcf}

## Our call have 19900 variants, 
## have format and sample filed
mutect2_vcf=${mutect2_dir}/${tumour_id}.final.vcf.gz
mutect2_variants=$(bcftools view -H ${mutect2_vcf} | wc -l)
echo "Mutect2 variants: $mutect2_variants"
# bcftools query -l ${mutect2_vcf}
# zgrep "^#CHROM" ${mutect2_vcf}

## Truth set, indel=1625, snv=39447
truth_vcf_indel=/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395/raw/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz
truth_indel_variants=$(bcftools view -H ${truth_vcf_indel} | wc -l)
echo "Truth indel variants: $truth_indel_variants"

truth_vcf_snv=/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395/raw/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz
truth_snv_variants=$(bcftools view -H ${truth_vcf_snv} | wc -l)
echo "Truth snv variants: $truth_snv_variants"

truth_bed=/home/zhonggr/projects/250224_DFSP_WES/data/benchmark/HCC1395/raw/High-Confidence_Regions_v1.2.bed

## Combine truth indel and snv into one VCF
export HGREF="${reference}"
conda activate hap

# hap_dir=${work_dir}/hap/mutect2_vs_truthmutect2
# mkdir -p ${hap_dir}

# hap.py \
#     ${truth_mutect2_vcf} \
#     ${mutect2_vcf} \
#     -r ${reference} \
#     -f ${truth_bed} \
#     -o ${hap_dir}/hap \
#     --verbose

hap_dir=${work_dir}/hap/mutect2_vs_truth_snv
mkdir -p ${hap_dir}


hap.py \
    ${truth_vcf_snv} \
    ${mutect2_vcf} \
    -r ${reference} \
    -f ${truth_bed} \
    -o ${hap_dir}/hap \
    --verbose
