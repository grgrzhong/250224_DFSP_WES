#!/bin/bash
#SBATCH --job-name=CNV_FACETS
#SBATCH --partition=amd
#SBATCH --time=4:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate facets

## Function to run CNV FACETS for one sample
run_cnv_facets() {
    local tumour_id="$1"
    local ref_dir="$2"
    local bam_dir="$3"
    local out_dir="$4"
    local dbsnp="$5"
    local dbsnp_index="$6"
    local defined_normal="$7"

    # Extract case_id from tumour_id (everything before the second hyphen)
    case_id=$(echo $tumour_id | cut -d'-' -f1,2)
    normal_id=${case_id}-N
    
    #Check if presence of paired normal samples
    if [ -d ${bam_dir}/${normal_id} ]; then

        meta_id=${tumour_id}_vs_${normal_id}

        echo "Processing : ${meta_id}"

        # Create temp directory for this sample
        mkdir -p "${out_dir}/${tumour_id}"
        
        prefix="${out_dir}/${tumour_id}/${tumour_id}"

        # Run FACETS with paired normal sample
        cnv_facets.R \
            --snp-tumour $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            --snp-normal $bam_dir/${normal_id}/${normal_id}_recalibrated.bam \
            --snp-vcf ${dbsnp} \
            --snp-nprocs 4 \
            --out ${prefix} \
            >& ${prefix}.facets.log
        
    else

        echo "No paired normal sample found, using manually definied normal sample ${defined_normal}"
        
        meta_id=${tumour_id}_vs_${defined_normal}
        
        # Create temp directory for this sample
        mkdir -p "${out_dir}/${tumour_id}"
        
        prefix="${out_dir}/${tumour_id}/${tumour_id}"
        
        # Run FACETS with manually definied normal sample
        cnv_facets.R \
            --snp-tumour $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            --snp-normal $bam_dir/${defined_normal}/${defined_normal}_recalibrated.bam \
            --snp-vcf ${dbsnp} \
            --snp-nprocs 1 \
            --out ${prefix} \
            >& ${prefix}.facets.log

    fi
}

export -f run_cnv_facets

## input parameters
ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
out_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets_test_normal"
out_dir="/tmp"
mkdir -p "${out_dir}"

dbsnp="${ref_dir}/dbSNP.vcf.gz"
dbsnp_index="${ref_dir}/dbSNP.vcf.gz.tbi"

## Define a default normal sample to use when paired normal is not available
defined_normal="DFSP-336-N"

## sample list
# tumour_samples=/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/sample_list.txt
tumour_samples=/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/tumour_only_samples.txt
# tumour_samples=/home/zhonggr/projects/250224_DFSP_WES/data/wes/sample_info/tumour_test1.txt

# Count the number of samples
num_sample=$(cat "$tumour_samples" | wc -l)
echo "Number of samples = $num_sample"

# Calculate number of parallel jobs
if [ $num_sample -ge 15 ]; then
    jobs=15
else
    jobs=$num_sample
fi

echo "Parallel jobs = $jobs"

## Run FACETS for each tumour sample
cat ${tumour_samples} | parallel \
    --jobs $jobs \
    --progress \
    run_cnv_facets {} "$ref_dir" "$bam_dir" "$out_dir" "$dbsnp" "$dbsnp_index" "$defined_normal"