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

# singularity shell \
#     --bind /home/zhonggr/projects/250224_DFSP_WES/data/:/data \
#     --bind /home/zhonggr/projects/250224_DFSP_WES/data/reference:/reference \
#     /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/cnv_facets.sif

# conda create -n facets 
# conda install cnv_facets

#

# bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
# data_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
# ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
# tumour_bam="${data_dir}/DFSP-336-T/DFSP-336-T_recalibrated.bam"
# tumour_bai="${data_dir}/DFSP-336-T/DFSP-336-T_recalibrated.bai"
# normal_bam="${data_dir}/DFSP-336-N/DFSP-336-N_recalibrated.bam"
# normal_bai="${data_dir}/DFSP-336-N/DFSP-336-N_recalibrated.bai"
# prefix="/home/zhonggr/projects/250224_DFSP_WES/data/test/cnv/facets/DFSP-336-T/DFSP-336-T"

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate facets

## Function to run CNV FACETS
cnv_facets_call() {
    local tumour_id="$1"
    local ref_dir="$2"
    local bam_dir="$3"
    local tmp_dir="$4"
    local out_dir="$5"
    local dbsnp="$6"
    local dbsnp_index="$7"
    local defined_normal="$8"

    # Extract case_id from tumour_id (everything before the second hyphen)
    case_id=$(echo $tumour_id | cut -d'-' -f1,2)
    normal_id=${case_id}-N
    
    
    #Check if presence of paired normal samples
    if [ -d ${bam_dir}/${normal_id} ]; then

        # meta_id=${tumour_id}_vs_${normal_id}
        meta_id=${tumour_id}

        echo "Processing sample: ${meta_id}"

        # Create temp directory for this sample
        mkdir -p "${tmp_dir}/${tumour_id}"
        
        prefix="${tmp_dir}/${tumour_id}/${meta_id}"

        # Run FACETS with paired normal sample
        cnv_facets.R \
            --snp-tumour $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            --snp-normal $bam_dir/${normal_id}/${normal_id}_recalibrated.bam \
            --snp-vcf ${dbsnp} \
            --snp-nprocs 4 \
            --out ${prefix} \
            >& ${prefix}.facets.log

        # Create output directory 
        cp ${tmp_dir}/${tumour_id} ${out_dir}/
        
    else

        echo "No paired normal sample found, using manually definied normal sample ${defined_normal}"
        # meta_id=${tumour_id}_vs_${defined_normal}
        meta_id=${tumour_id}
        
        # Create temp directory for this sample
        mkdir -p "${tmp_dir}/${tumour_id}"
        
        prefix="${tmp_dir}/${tumour_id}/${meta_id}"
        
        # Run FACETS with manually definied normal sample
        cnv_facets.R \
            --snp-tumour $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            --snp-normal $bam_dir/${defined_normal}/${defined_normal}_recalibrated.bam \
            --snp-vcf ${dbsnp} \
            --snp-nprocs 4 \
            --out ${prefix} \
            >& ${prefix}.facets.log
        
        mv ${tmp_dir}/${tumour_id} ${out_dir}/
    fi
}

export -f cnv_facets_call

## Set directories
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
out_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets_test_normal"
mkdir -p "${out_dir}"
dbsnp="${ref_dir}/dbSNP.vcf.gz"
dbsnp_index="${ref_dir}/dbSNP.vcf.gz.tbi"

# Create temporary directory on local filesystem
tmp_dir="/tmp/facets"
mkdir -p "${tmp_dir}"

# Define a default normal sample to use when paired normal is not available
defined_normal="DFSP-336-N"  # Set your default normal sample here

# sample list
# tumour_samples=/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/sample_list.txt
tumour_samples=/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/tumour_only_samples.txt

## Run FACETS for each tumour sample
parallel_jobs=10
cat ${tumour_samples} | parallel \
    --jobs $parallel_jobs \
    --progress \
    cnv_facets_call {} "$ref_dir" "$bam_dir" "$tmp_dir" "$out_dir" "$dbsnp" "$dbsnp_index" "$defined_normal"