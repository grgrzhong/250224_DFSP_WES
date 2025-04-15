#!/bin/bash

# Setup the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NFX_OPTS="-Xms=512m -Xmx=8g"
export NXF_ANSI_LOG=true

export NXF_LOG_FILE="${PWD}/test/.nextflow.log"
rm -f ${NXF_LOG_FILE}

# Test the mutation_calling workflow
nextflow run modules/variant_calling/vcf2maf/vcf2maf/test/main.nf.test \
    -profile local \
    -work-dir test/work