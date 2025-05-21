#!/bin

# singularity shell \
#     --bind /home/zhonggr/projects/250224_DFSP_WES/data/:/data \
#     --bind /home/zhonggr/projects/250224_DFSP_WES/data/reference:/reference \
#     /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/cnv_facets.sif

data_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
tumour_bam="${data_dir}/DFSP-336-T/DFSP-336-T_recalibrated.bam"
tumour_bai="${data_dir}/DFSP-336-T/DFSP-336-T_recalibrated.bai"
normal_bam="${data_dir}/DFSP-336-N/DFSP-336-N_recalibrated.bam"
normal_bai="${data_dir}/DFSP-336-N/DFSP-336-N_recalibrated.bai"
dbsnp="${ref_dir}/dbSNP.vcf.gz"
dbsnp_index="${ref_dir}/dbSNP.vcf.gz.tbi"

# Create temporary directory on local filesystem
temp_dir="/tmp/facets_temp_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${temp_dir}"
prefix="${temp_dir}/DFSP-336-T"

# prefix="/home/zhonggr/projects/250224_DFSP_WES/data/test/cnv/facets/DFSP-336-T/DFSP-336-T"

# Run FACETS
cnv_facets.R \
    --snp-normal ${normal_bam} \
    --snp-tumour ${tumour_bam} \
    --snp-vcf ${dbsnp} \
    --snp-nprocs 4 \
    --out ${prefix}

# Copy results to final destination if successful
if [ $? -eq 0 ]; then
    final_dir="/home/zhonggr/projects/250224_DFSP_WES/data/test/cnv/facets/DFSP-336-T"
    mkdir -p "${final_dir}"
    cp "${temp_dir}"/* "${final_dir}/"
    echo "Results copied to ${final_dir}"
fi