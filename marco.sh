#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Setup the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS="-Xms512m -Xmx8g"
# export NXF_LOG_FILE="${PWD}/.nextflow.log"
# Run marco
export NXF_LOG_FILE="${PWD}/test/marco.nextflow.log"
# rm -rf ${NXF_LOG_FILE}
nextflow run subworkflows/mutation_calling/mutect2_call_only.nf \
    -profile hpc \
    --input /lustre1/g/path_my/samplesheet.csv \
    --panel_of_normals /lustre1/g/path_my/30X_FT/SNVs/Mutect2-Call/PON-Mutect/pon.vcf.gz \
    --panel_of_normals_tbi /lustre1/g/path_my/30X_FT/SNVs/Mutect2-Call/PON-Mutect/pon.vcf.gz.tbi \
    --assay_type wgs \
    -work-dir test/work \
    --outdir test/results

# export NXF_LOG_FILE="${PWD}/marco2.nextflow.log"
# rm -f ${NXF_LOG_FILE}
# nextflow run workflows/mutation_calling/main.nf \
#     -work-dir work \
#     -profile hpc \
#     -resume \
#     --input /lustre1/g/path_my/samplesheet_P.csv \
#     --panel_of_normals /lustre1/g/path_my/30X_cfdna/SNVs/Mutect2-Call/PON-Mutect/pon.vcf.gz \
#     --step mutect2 \
#     --outdir test/macro2