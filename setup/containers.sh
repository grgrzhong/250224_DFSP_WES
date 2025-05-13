#!/bin/bash

# Download containers from the Seqera community registry
# using Singularity. The containers are used for the CNV_FACETS workflow.
# https://seqera.io/containers/

container_dir="/home/zhonggr/projects/250224_DFSP_WES/containers/singularity"
mkdir -p ${container_dir}

singularity pull --force --dir ${container_dir} bamtools.sif oras://community.wave.seqera.io/library/bamtools:2.5.2--ec8f9631801f9901
singularity pull --force --dir ${container_dir} bcftools.sif oras://community.wave.seqera.io/library/bcftools:1.21--21573c18b3ab6bcb
singularity pull --force --dir ${container_dir} bwa.sif oras://community.wave.seqera.io/library/bwa_samtools:3f723aba7e77ee82
singularity pull --force --dir ${container_dir} cnv_facets.sif oras://community.wave.seqera.io/library/cnv_facets:0.16.1--279b2b94c5e037b9
singularity pull --force --dir ${container_dir} cnvkit.sif oras://community.wave.seqera.io/library/cnvkit:0.9.12--8f4ba584e385f393
singularity pull --force --dir ${container_dir} controlfreec.sif oras://community.wave.seqera.io/library/control-freec:11.6--b8904bfc98b3c9ba
singularity pull --force --dir ${container_dir} delly.sif oras://community.wave.seqera.io/library/delly:1.3.3--6a757953473b323c
singularity pull --force --dir ${container_dir} ensemblvep.sif oras://community.wave.seqera.io/library/ensembl-vep:113.4--ca721b5ce97a9cdc
singularity pull --force --dir ${container_dir} fastp.sif oras://community.wave.seqera.io/library/fastp:0.24.0--0397de619771c7ae
singularity pull --force --dir ${container_dir} fastqc.sif oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960
singularity pull --force --dir ${container_dir} freebayes.sif oras://community.wave.seqera.io/library/freebayes:1.3.9--248f50c2ec249a6f
singularity pull --force --dir ${container_dir} gatk4.sif oras://community.wave.seqera.io/library/gatk4:4.5.0.0--7e714e2ee9c2e11e
singularity pull --force --dir ${container_dir} manta.sif oras://community.wave.seqera.io/library/manta:1.6.0--765006f7a2a42f9a
singularity pull --force --dir ${container_dir} multiqc.sif oras://community.wave.seqera.io/library/multiqc:1.28--d466e41d58d6d704
singularity pull --force --dir ${container_dir} pysam.sif oras://community.wave.seqera.io/library/pysam_python:65ddcfde3934cade
singularity pull --force --dir ${container_dir} samtools.sif oras://community.wave.seqera.io/library/samtools:1.21--84c9d77c3901e90b
singularity pull --force --dir ${container_dir} tabix.sif oras://community.wave.seqera.io/library/tabix:1.11--dba91ce963b95ef9
singularity pull --force --dir ${container_dir} vcf2maf.sif oras://community.wave.seqera.io/library/ensembl-vep_vcf2maf:df9fb21816c1659e



