#!/bin/bash

# Set the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NFX_OPTS="-Xms=512m -Xmx=8g"
export NXF_ANSI_LOG=false

# Set the Nextflow log file and remove the existing log file
export NXF_LOG_FILE="${PWD}/logs/nextflow/nextflow.log"

rm -f ${NXF_LOG_FILE}

# Run the Nextflow workflow
nextflow run workflows/mutation_calling/main.nf \
    --input /home/zhonggr/projects/250224_DFSP_WES/data/WES/csv/test.csv \
    --outdir ./data/WES
