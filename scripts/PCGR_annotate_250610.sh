#!/bin/bash
#SBATCH --job-name=PCGR_Annotation
#SBATCH --partition=amd
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

## =============================================================================
## PCGR Annotation (https://sigven.github.io/pcgr/index.html)
## Authors: Zhong Guorui
## Date: 2025-06-10
## Description: This script was adapted from "M:\Scripts\DNA analysis\PCGR_annotate.sh"
## Added prallelization and other features for PCGR annotation
## =============================================================================

## Project settings
export project_dir="/mnt/f/projects/250224_DFSP_WES"
export module_dir="${project_dir}/bin"
export bam_dir="${project_dir}/data/wes/bam"
export mutect2_dir="${project_dir}/data/wes/mutect2"

## Reference settings
export ref_dir="/mnt/m/Reference"
export ref_data_dir="${ref_dir}/PCGR_reference/20250314"
export vep_dir="${ref_dir}/VEP_cache"
export panel_of_normals="${ref_dir}/WES/DFSP/PON-Mutect/pon.vcf.gz"

## Download PCGR singularity container if not already present
container_dir="${project_dir}/containers"
singularity pull --force --dir ${container_dir} pcgr-2.2.1.sif oras://ghcr.io/sigven/pcgr:2.2.1.singularity

## Output directory for PCGR results
work_dir="${project_dir}/data/wes/PCGR"

## Sample list to run PCGR annotation
sample_list=$(ls $mutect2_dir)

## Parallel jobs
num_sample=$(echo "$sample_list" | wc -l)
if [ "${num_sample}" -ge 15 ]; then    
    jobs=15
else
    jobs=$num_sample
fi

# tumour_id="DFSP-001-T"
tumour_id="DFSP-028-T"
# tumour_id="DFSP-031-T"
## Function to run PCGR annotation
pcgr_annotation() {
    local tumour_id=$1
    local mutect2_dir=$2
    local bam_dir=$3
    local ref_data_dir=$4
    local vep_dir=$5
    local work_dir=$6
    local panel_of_normals=$7
    local container_dir=$8

    ## Find the matched normal sample
    patient_id=$(echo "${tumour_id}" | cut -d'-' -f1,2)
    normal_id="${patient_id}-N"

    echo "$(date +"%F") $(date +"%T")" "Processing sample: ${tumour_id} for patient: ${patient_id}"

    ## Check if normal sample exists and set is_paired flag
    if [ ! -d "${bam_dir}/${normal_id}" ]; then
        is_paired="false"
    else
        is_paired="true"
    fi

    ## Check if input VCF exists
    input_vcf="${mutect2_dir}/${tumour_id}/${tumour_id}.final.vcf.gz"

    if [ ! -f "${input_vcf}" ]; then
        echo "Input VCF not found: ${input_vcf}"
        return 1
    fi
    
    ## Create output directory for this sample
    output_dir="${work_dir}/${tumour_id}"
    mkdir -p "${output_dir}"

    ## Define reformatted VCF path
    reformatted_vcf="${output_dir}/${tumour_id}.reformatted.vcf.gz"

    ## ========================================================================
    ## Tumour-Normal PCGR Annotation
    ## ========================================================================
    if [ "${is_paired}" = "true" ]; then

        echo "$(date +"%F") $(date +"%T")" "Normal sample not found for ${tumour_id}. Running tumor-only analysis ..."

        ## Reformat VCF 
        singularity exec \
            --bind "${mutect2_dir}:${mutect2_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind "/tmp:/tmp" \
            "${container_dir}/pysam-0.23.2.sif" \
            python "${module_dir}/pcgr_reformat_vcf_tumour_normal.py" \
                --input "${input_vcf}" \
                --output "${reformatted_vcf}"
        
        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            "${container_dir}/tabix-1.11.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${ref_data_dir}:${ref_data_dir}" \
            --bind "${vep_dir}:${vep_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${work_dir}:${work_dir}" \
            "${container_dir}/pcgr-2.2.1.sif" \
            pcgr \
            --input_vcf "${reformatted_vcf}" \
            --vep_dir "${vep_dir}" \
            --refdata_dir "${ref_data_dir}" \
            --output_dir "${output_dir}" \
            --genome_assembly grch38 \
            --sample_id "${tumour_id}" \
            --assay WES \
            --effective_target_size_mb 34 \
            --tumor_dp_tag TDP \
            --tumor_af_tag TAF \
            --control_dp_tag NDP \
            --control_af_tag NAF \
            --tumor_dp_min 20 \
            --tumor_af_min 0.05 \
            --control_dp_min 10 \
            --control_af_max 0.01 \
            --estimate_tmb \
            --tmb_dp_min 20 \
            --tmb_af_min 0.05 \
            --estimate_msi \
            --estimate_signatures \
            --vcf2maf \
            --ignore_noncoding \
            --force_overwrite \
            >& "${output_dir}/pcgr.log"

    ## ========================================================================
    ## Tumour-Only PCGR Annotation
    ## ========================================================================
    else
        
        echo "$(date +"%F") $(date +"%T")" "Normal sample not found for ${tumour_id}. Running tumor-only analysis ..."
        
        ## Reformat VCF 
        singularity exec \
            --bind "${mutect2_dir}:${mutect2_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${module_dir}:${module_dir}" \
            --bind "/tmp:/tmp" \
            "${container_dir}/pysam-0.23.2.sif" \
            python "${module_dir}/pcgr_reformat_vcf_tumour_only.py" \
                --input "${input_vcf}" \
                --output "${reformatted_vcf}"
        
        ## Index it
        rm -f "${reformatted_vcf}.tbi"
        singularity exec \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            "${container_dir}/tabix-1.11.sif" \
            tabix -p vcf "${reformatted_vcf}"
        
        ## Run PCGR annotation
        singularity exec \
            --bind "${ref_data_dir}:${ref_data_dir}" \
            --bind "${vep_dir}:${vep_dir}" \
            --bind "${bam_dir}:${bam_dir}" \
            --bind "${mutect2_dir}:${mutect2_dir}" \
            --bind "${work_dir}:${work_dir}" \
            --bind "${output_dir}:${output_dir}" \
            --bind "${module_dir}:${module_dir}" \
            "${container_dir}/pcgr-2.2.1.sif" \
            pcgr \
            --input_vcf "${reformatted_vcf}" \
            --vep_dir "${vep_dir}" \
            --refdata_dir "${ref_data_dir}" \
            --output_dir "$output_dir" \
            --genome_assembly grch38 \
            --sample_id "${tumour_id}" \
            --assay WES \
            --effective_target_size_mb 34 \
            --tumor_only \
            --tumor_dp_tag TDP \
            --tumor_af_tag TAF \
            --tumor_dp_min 20 \
            --tumor_af_min 0.05 \
            --estimate_tmb \
            --tmb_dp_min 20 \
            --tmb_af_min 0.05 \
            --estimate_msi \
            --estimate_signatures \
            --vcf2maf \
            --ignore_noncoding \
            --force_overwrite \
            >& "${output_dir}/pcgr_annotation.log"

    fi
}


