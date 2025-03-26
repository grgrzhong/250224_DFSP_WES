#!/bin/bash

#############################################################################
# Project configuration for basic path variables
# Authors: Zhong Guorui
# Date: 2024-03-12
#############################################################################

# Prjoect and working directory
export BASE_DIR="/lustre1/g/path_my"
export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export PROJECT_DIR="250224_DFSP_WES"
# export WORK_DIR="data/WES/DFSP"
export WORK_DIR="data/WES/SARC"

export SAMPLE_SHEET=$BASE_DIR/$PROJECT_DIR/$WORK_DIR/samplesheet.csv

# Extract sample information
export CASE_COUNT=$(tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f1 | sort -u | wc -l)
export TUMOUR_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f2 | grep -E "\-T")
export NORMAL_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | cut -d',' -f2 | grep -E "\-N")
export TUMOUR_COUNT=$(echo "$TUMOUR_SAMPLES" | wc -w)
export NORMAL_COUNT=$(echo "$NORMAL_SAMPLES" | wc -w)

# Setup PARALLEL_JOBS based on case count
if [ "$CASE_COUNT" -gt 30 ]; then
    export PARALLEL_JOBS=30
else
    export PARALLEL_JOBS=$CASE_COUNT
fi

# NXF_OPTS='-Xms1g -Xmx4g'

# Setup memory usage for GATK java options
# TOTAL_MEM=$(free -g | awk '/^Mem:/{print $2}')
# export AVAIL_MEM=$(($TOTAL_MEM * 80 / 100 / $PARALLEL_JOBS))

# if [ "$AVAIL_MEM" -lt 4 ]; then
#     export AVAIL_MEM=4 # Set minimum memory threshold
# fi

# export MAX_MEM=$(($TOTAL_MEM * 8 / 10 /))

# Load general utility functions
source "${BASE_DIR}/${PROJECT_DIR}/modules/utils/utils.sh"
