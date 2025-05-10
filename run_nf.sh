#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Setup the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS="-Xms512m -Xmx8g"
# export NXF_LOG_FILE="${PWD}/.nextflow.log"
rm -f .nextflow.log*

nextflow run subworkflows/mutation_calling/mutect2_call.nf \
    -profile hpc \
    --input /home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test1.csv

## Test the cnv_facets
# nextflow run subworkflows/mutation_calling/cnv_facets_old.nf \
#     -profile hpc \
#     -resume \
#     --outdir data/wes

# Test the mutect2
# export NXF_LOG_FILE="${PWD}/.nextflow.log"
# rm -f ${NXF_LOG_FILE}
# nextflow run subworkflows/mutation_calling/mutect2_call.nf \
#     -profile hpc \
#     --input /home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test1.csv \
#     --outdir test/mutect2

# # Run the workflow local
# nextflow run workflows/mutation_calling/main.nf \
#     -profile local \
#     -resume \
#     --input /home/zhonggr/projects/250224_DFSP_WES/data/sarc/csv/samplesheet.csv \
#     --outdir data/sarc

# run testing
# rm -f ${NXF_LOG_FILE}
# nextflow run workflows/mutation_calling/main.nf \
#     -profile hpc \
#     -resume \
#     -work-dir work \
#     --step mutect2 \
#     --outdir results

# Run marco
# export NXF_LOG_FILE="${PWD}/marco.nextflow.log"
# rm -f ${NXF_LOG_FILE}
# nextflow run workflows/mutation_calling/main.nf \
#     -work-dir work \
#     -profile hpc \
#     -resume \
#     --input /lustre1/g/path_my/samplesheet.csv \
#     --panel_of_normals /lustre1/g/path_my/30X_FT/SNVs/Mutect2-Call/PON-Mutect/pon.vcf.gz \
#     --step mutect2 \
#     --outdir test/macro

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