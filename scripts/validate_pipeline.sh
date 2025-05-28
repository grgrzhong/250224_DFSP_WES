#!/bin/bash
#SBATCH --job-name=PreProcess_WES
#SBATCH --partition=amd
#SBATCH --time=96:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall


# Define reference files
export ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
export reference="${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed"
export dbsnp="${ref_dir}/Population_database/dbSNP.vcf.gz"
export bait_interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list"
export target_interval="${ref_dir}/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list"



sample=NA12878
fastq_1=/home/zhonggr/projects/250224_DFSP_WES/data/reference/benchmark/NA12878/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
fastq_2=/home/zhonggr/projects/250224_DFSP_WES/data/reference/benchmark/NA12878/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# zcat ${fastq_1} | head -n 15
# zcat ${fastq_2} | head -n 15

work_dir=/home/zhonggr/projects/250224_DFSP_WES/data/reference/benchmark/NA12878
    
# Check if input files exist
if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
    echo "ERROR: Input FASTQ files not found for sample $sample" >&2
    return 1
fi

# Create output directories
fastq_trimmed_dir="${work_dir}/preprocessing/fastq_trimmed"
mkdir -p ${fastq_trimmed_dir}/${sample}

bam_dir="${work_dir}/preprocessing/bam"
mkdir -p ${bam_dir}/${sample}

echo "Processing sample:        $sample"
echo "FASTQ1:                   $fastq_1"
echo "FASTQ2:                   $fastq_2"
echo "Working directory:        $work_dir"
echo "BAM directory:            $bam_dir"
echo "Fastq-trimmed directory:  $bam_dir"


# Experimental UMI options
umi_opts=""
use_umi=false
umi_loc="per_read"
umi_len=8
if [[ "${use_umi}" == "true" ]]; then
    umi_opts="--umi --umi_loc ${umi_loc} --umi_len ${umi_len}"
fi


# adapter_1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
# adapter_2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

adapter_1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
adapter_2="CTGTCTCTTATACACATCTGACGCTGCCGACGA"

# Step 1: Fastp trimming
echo $(date +"%F") $(date +"%T") "Running Fastp trimming for $sample..."
fastp \
    -i ${fastq_1} \
    -I ${fastq_2} \
    -o ${fastq_trimmed_dir}/${sample}/${sample}_trimmed_1.fastq.gz \
    -O ${fastq_trimmed_dir}/${sample}/${sample}_trimmed_2.fastq.gz \
    --detect_adapter_for_pe \
    --adapter_sequence=${adapter_1} \
    --adapter_sequence_r2=${adapter_2} \
    --trim_poly_g \
    --trim_poly_x \
    ${umi_opts} \
    -j ${fastq_trimmed_dir}/${sample}/${sample}.json \
    -h ${fastq_trimmed_dir}/${sample}/${sample}.html \
    -w 8

if [[ $? -ne 0 ]]; then
    echo "ERROR: Fastp trimming failed for ${sample}" >&2
    return 1
fi

# Step 2: BWA Alignment
echo $(date +"%F") $(date +"%T") "Aligning ${sample} to reference genome..."

kit_name="Nextera_Rapid_Capture_Exome"
platform="ILLUMINA"
platform_model="HISEQ2500" 

bwa mem -M -t 16 \
    -R "@RG\tID:${sample}\tLB:${kit_name}\tPL:${platform}\tPM:${platform_model}\tSM:${sample}\tPU:NA" \
    ${reference} \
    ${fastq_trimmed_dir}/${sample}/${sample}_trimmed_1.fastq.gz \
    ${fastq_trimmed_dir}/${sample}/${sample}_trimmed_2.fastq.gz > ${bam_dir}/${sample}/${sample}.sam

if [[ $? -ne 0 ]]; then
    echo "ERROR: BWA alignment failed for ${sample}" >&2
    return 1
fi

samtools view -Sb ${bam_dir}/${sample}/${sample}.sam > ${bam_dir}/${sample}/${sample}.bwa.bam


# # Step 3: Extract and tag UMI, remove as 
# echo $(date +"%F") $(date +"%T") "Extract and tag UMI for ${sample}..."
# python /home/zhonggr/projects/250224_DFSP_WES/bin/python/tag_umi.py \
#     -I ${bam_dir}/${sample}/${sample}.bwa.bam \
#     -O ${bam_dir}/${sample}/${sample}.umi.bam

if [[ $? -ne 0 ]]; then
    echo "ERROR: UMI tagging failed for ${sample}" >&2
    return 1
fi

# Step 4: Sort by coordinate
echo $(date +"%F") $(date +"%T") "Sorting BAM file for ${sample}..."
samtools sort \
    ${bam_dir}/${sample}/${sample}.bwa.bam \
    -o ${bam_dir}/${sample}/${sample}.sorted.bam

if [[ $? -ne 0 ]]; then
    echo "ERROR: BAM sorting failed for ${sample}" >&2
    return 1
fi

# Step 5: Mark duplicates
echo $(date +"%F") $(date +"%T") "Marking duplicates for ${sample}..."
gatk --java-options -Xmx4g MarkDuplicates \
    -I ${bam_dir}/${sample}/${sample}.sorted.bam \
    -M ${bam_dir}/${sample}/${sample}.metrics.txt \
    -O ${bam_dir}/${sample}/${sample}.marked.bam \
    # Remove: --BARCODE_TAG "RX"  # No UMIs available

if [[ $? -ne 0 ]]; then
    echo "ERROR: Marking duplicates failed for ${sample}" >&2
    return 1
fi

# Step 6: Index BAM
echo $(date +"%F") $(date +"%T") "Indexing BAM file for ${sample}..."
samtools index ${bam_dir}/${sample}/${sample}.marked.bam

if [[ $? -ne 0 ]]; then
    echo "ERROR: BAM indexing failed for ${sample}" >&2
    return 1
fi

# Step 7: Base recalibration
echo $(date +"%F") $(date +"%T") "Running base recalibration for ${sample}..."
gatk BaseRecalibrator \
    -I ${bam_dir}/${sample}/${sample}.marked.bam \
    -R ${reference} \
    -L ${interval} \
    -O ${bam_dir}/${sample}/${sample}.recal.table \
    --known-sites ${dbsnp}

if [[ $? -ne 0 ]]; then
    echo "ERROR: Base recalibration failed for ${sample}" >&2
    return 1
fi

# Step 8: Apply BQSR
echo $(date +"%F") $(date +"%T") "Applying BQSR for ${sample}..."
gatk ApplyBQSR \
    -I ${bam_dir}/${sample}/${sample}.marked.bam \
    -O ${bam_dir}/${sample}/${sample}.recal.bam \
    -L ${interval} \
    -bqsr ${bam_dir}/${sample}/${sample}.recal.table \
    --create-output-bam-md5

if [[ $? -ne 0 ]]; then
    echo "ERROR: Applying BQSR failed for ${sample}" >&2
    return 1
fi

# Step 9: Collect metrics
echo $(date +"%F") $(date +"%T") "Collecting HsMetrics for ${sample}..."
gatk CollectHsMetrics \
    -I ${bam_dir}/${sample}/${sample}.recal.bam \
    -O ${bam_dir}/${sample}/${sample}.hsmetrics.txt \
    -R ${reference} \
    -BI ${bait_interval} \
    -TI ${target_interval}

if [[ $? -ne 0 ]]; then
    echo "ERROR: Collecting HsMetrics failed for ${sample}" >&2
    return 1
fi

# Step 10: Generate alignment stats
echo $(date +"%F") $(date +"%T") "Generating alignment stats for ${sample}..."
bamtools stats \
    -in ${bam_dir}/${sample}/${sample}.recal.bam \
    > ${bam_dir}/${sample}/${sample}.alnstat.txt

if [[ $? -ne 0 ]]; then
    echo "ERROR: Generating alignment stats failed for ${sample}" >&2
    return 1
fi

# Clean up intermediate files
echo $(date +"%F") $(date +"%T") "Cleaning up intermediate files for ${sample}..."
rm -f ${cur_dir}/${sample}.bwa.bam \
    ${cur_dir}/${sample}.umi.bam \
    ${cur_dir}/${sample}.sorted.bam \
    ${cur_dir}/${sample}.marked.bam \
    ${cur_dir}/${sample}.marked.bam.bai

echo $(date +"%F") $(date +"%T") "Completed processing sample: ${sample}"


# Create output directories
# export fastq_trimmed_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/fastq_trimmed"
# export bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
export fastq_trimmed_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/test/fastq_trimmed"
export bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/test/recalibrated"

mkdir -p ${fastq_trimmed_dir}
mkdir -p ${bam_dir}

# CSV file with sample information
# export input_csv="/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/fastq_all_samples.csv"
export input_csv="/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/fastq_test1.csv"

# Create a temporary queue file to handle special characters in paths
temp_file=$(mktemp)

# The CSV has format: "patient","sample","status","fastq_1","fastq_2"
tail -n +2 ${input_csv} | while IFS=, read -r patient sample status fastq_1 fastq_2; do
    # Remove quotes if present
    sample=$(echo "$sample" | tr -d '"')
    fastq_1=$(echo "$fastq_1" | tr -d '"')
    fastq_2=$(echo "$fastq_2" | tr -d '"')
    
    echo "${sample}|${fastq_1}|${fastq_2}" >> "${temp_file}"
done

# Count samples
num_samples=$(wc -l < "${temp_file}")

# Calculate number of parallel jobs
if [ ${num_samples} -ge 15 ]; then
    jobs=15
else
    jobs=${num_samples}
fi

## Print configuration
echo "=========================================================================="
echo "Preprocessing configuration:"
echo "=========================================================================="
echo "Input CSV:             ${input_csv}"
echo "Reference directory:   ${ref_dir}"
echo "Reference:             ${reference}"
echo "Interval:              ${interval}"
echo "dbSNP:                 ${dbsnp}"
echo "Bait interval:         ${bait_interval}"
echo "Target interval:       ${target_interval}"
echo "Trimmed fastq dir:     ${fastq_trimmed_dir}"
echo "BAM directory:         ${bam_dir}"
echo "Number of samples      ${num_samples}"
echo "Parallel jobs          ${jobs}"

# Process with parallel
cat "${temp_file}" | parallel \
    --jobs ${jobs} \
    --progress \
    --colsep '|' \
    run_preprocessing {1} {2} {3} "$fastq_trimmed_dir" "$bam_dir" "$reference" "$interval" "$dbsnp" "$bait_interval" "$target_interval"

# Clean up
rm -f "${temp_file}"

# Run FastQC on trimmed fastq files
# echo $(date +"%F") $(date +"%T") "Running FastQC on trimmed files..."
# mkdir -p ${fastq_trimmed_dir}/FastQC
# fastqc ${fastq_trimmed_dir}/*.fastq.gz \
#     -o ${fastq_trimmed_dir}/FastQC \
#     --memory 8192 -t 8

echo "All samples processed successfully."