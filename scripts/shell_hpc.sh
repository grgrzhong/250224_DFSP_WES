#!/bin/bash

#############################################################################
# Interactive HPC Session Request Script
# Description: Request resources for interactive usage on HPC
# Version: 1.0.0
#############################################################################

# Default settings
CPUS=16
MEM=64
TIME="10:00:00"
PARTITION="amd"
QOS="normal"
X11=false
JUPYTER=false
JUPYTER_PORT=8888

# Using srun for interactive session
# srun --partition=${PARTITION} \
#     --qos=${QOS} \
#     --nodes=1 \
#     --ntasks=1 \
#     --cpus-per-task=${CPUS} \
#     --mem=${MEM}G \
#     --time=${TIME} \
#     --pty bash -i

# Using salloc for interactive session
salloc --partition=${PARTITION} \
      --qos=${QOS} \
      --nodes=1 \
      --ntasks=1 \
      --cpus-per-task=${CPUS} \
      --mem=${MEM}G \
      --time=${TIME}

# Conditional X11 forwarding
if [ "$X11" = true ]; then
    ssh -X $SLURM_NODELIST
fi