#!/bin/bash

# Download containers from the Seqera community registry
# using Singularity. The containers are used for the CNV_FACETS workflow.
# https://seqera.io/containers/

container_dir="/home/zhonggr/projects/250224_DFSP_WES/containers"
mkdir -p ${container_dir}/singularity

singularity pull --dir ${container_dir}/singularity bamtools.sif oras://community.wave.seqera.io/library/bamtools:2.5.2--ec8f9631801f9901
singularity pull --dir ${container_dir}/singularity bcftools.sif oras://community.wave.seqera.io/library/bcftools:1.21--21573c18b3ab6bcb
singularity pull --dir ${container_dir}/singularity bwa.sif oras://community.wave.seqera.io/library/bwa_samtools:3f723aba7e77ee82
singularity pull --dir ${container_dir}/singularity cnv_facets.sif oras://community.wave.seqera.io/library/cnv_facets:0.16.1--279b2b94c5e037b9
singularity pull --dir ${container_dir}/singularity cnvkit.sif oras://community.wave.seqera.io/library/cnvkit:0.9.12--8f4ba584e385f393
singularity pull --dir ${container_dir}/singularity controlfreec.sif oras://community.wave.seqera.io/library/control-freec:11.6--b8904bfc98b3c9ba
singularity pull --dir ${container_dir}/singularity delly.sif oras://community.wave.seqera.io/library/delly:1.3.3--6a757953473b323c
singularity pull --dir ${container_dir}/singularity ensemblvep.sif oras://community.wave.seqera.io/library/ensembl-vep:113.4--ca721b5ce97a9cdc
singularity pull --dir ${container_dir}/singularity fastp.sif oras://community.wave.seqera.io/library/fastp:0.24.0--0397de619771c7ae
singularity pull --dir ${container_dir}/singularity fastqc.sif oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960
singularity pull --dir ${container_dir}/singularity freebayes.sif oras://community.wave.seqera.io/library/freebayes:1.3.9--248f50c2ec249a6f
singularity pull --dir ${container_dir}/singularity gatk4.sif oras://community.wave.seqera.io/library/gatk4:4.5.0.0--7e714e2ee9c2e11e
singularity pull --dir ${container_dir}/singularity manta.sif oras://community.wave.seqera.io/library/manta:1.6.0--765006f7a2a42f9a
singularity pull --dir ${container_dir}/singularity multiqc.sif oras://community.wave.seqera.io/library/multiqc:1.28--d466e41d58d6d704
singularity pull --dir ${container_dir}/singularity pysam.sif oras://community.wave.seqera.io/library/pysam_python:65ddcfde3934cade
singularity pull --dir ${container_dir}/singularity samtools.sif oras://community.wave.seqera.io/library/samtools:1.21--84c9d77c3901e90b
singularity pull --dir ${container_dir}/singularity vcf2maf.sif oras://community.wave.seqera.io/library/vcf2maf:1.6.22--478125c4e2c927b2



