#!/bin/bash

# Download containers from the Seqera community registry
# using Singularity. The containers are used for the CNV_FACETS workflow.
# https://seqera.io/containers/

container_dir="/home/zhonggr/projects/250224_DFSP_WES/containers"
mkdir -p ${container_dir}/singularity

singularity pull --dir ${container_dir}/singularity delly.sif oras://community.wave.seqera.io/library/delly:1.3.3--6a757953473b323c
singularity pull --dir ${container_dir}/singularity manta.sif oras://community.wave.seqera.io/library/manta:1.6.0--765006f7a2a42f9a
singularity pull --dir ${container_dir}/singularity cnv_facets.sif oras://community.wave.seqera.io/library/cnv_facets:0.16.1--279b2b94c5e037b9
