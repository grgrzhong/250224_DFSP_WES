#!/bin/bash

####
# Download the VEP cache file

vep_data_dir="$PWD/data/reference/ensembl_vep"
mkdir -p vep_data_dir
cd ${vep_data_dir}
wget -O http://ftp.ensembl.org/pub/release-113/variation/vep/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf *.tar.gz
rm homo_sapiens_vep_110_GRCh38.tar.gz # Optionally, remove the tar.gz file after extraction
echo "VEP cache downloaded and extracted to $VEP_DATA_DIR"


## ======================== Gene Annotation for GRCh38 ========================
gene_annotation_dir="$PWD/data/reference/gene_annotation"
mkdir -p ${gene_annotation_dir}
cd ${gene_annotation_dir}

# GENCODE 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
gunzip gencode.v43.basic.annotation.gtf.gz

# Ensembl
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz

## ======================== repeatmasker regions ========================
# url="https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1237006675_NR5drEgvZ85edZVWA5an0D0VQBC4&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_rmsk&hgta_ctDesc=table+browser+query+on+rmsk&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED"
# wget  -O- -q  "$url" | 
# grep -v "#" | 
# sort -k1,1 -k2,2n -k3,3n -t$'\t' | 
# bgzip -c > repeatmasker.GRCh38.bed.gz
# tabix -p bed repeatmasker.GRCh38.bed.gz

# Sort the repeatmasker bed file and create index
# gunzip ${ref_dir}/RepeatMasker.bed.gz
# bgzip ${ref_dir}/RepeatMasker.bed
# sort -k1,1 -k2,2n ${ref_dir}/RepeatMasker.bed -o ${ref_dir}/RepeatMasker.sorted.bed
# bgzip -c ${ref_dir}/RepeatMasker.sorted.bed > ${ref_dir}/RepeatMasker.sorted.bed.gz
# tabix -p bed ${ref_dir}/RepeatMasker.sorted.bed.gz

# # Sort the blacklist bed file and create index
# gunzip ${ref_dir}/blacklist.bed.gz
# bgzip ${ref_dir}/blacklist.bed
# sort -k1,1 -k2,2n ${ref_dir}/blacklist.bed -o ${ref_dir}/blacklist.sorted.bed
# bgzip -c ${ref_dir}/blacklist.sorted.bed > ${ref_dir}/blacklist.sorted.bed.gz
# tabix -p bed ${ref_dir}/blacklist.sorted.bed.gz

## InterVar
intervar_dir="$PWD/data/reference/InterVar-2.2.1"
wget -P ${intervar_dir}/intervardb https://www.omim.org/static/omim/data/mim2gene.txt

## annovar
annovar_dir="$PWD/data/annovar"
wget -P ${annovar_dir}/annovar http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -xzf ${annovar_dir}/annovar.latest.tar.gz -C ${annovar_dir}/annovar

wget https://repo.bioserver.ieo.it/dima/hg38/humandb/hg38_ALL.sites.2015_08.txt


# Download COSMIC (requires license, but a free subset is available)
cosmic_dir="$PWD/data/reference/cosmic"
wget https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/vXX/CosmicMutantExport.tsv.gz

# Download ClinVar (public)
clinvar_dir="$PWD/data/reference/clinvar"
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# Download OncoKB (free version available)
oncokb_dir="$PWD/data/reference/OncoKB"
mkdir -p ${oncokb_dir}
wget https://www.oncokb.org/api/v1/utils/allAnnotatedVariants.tsv -O oncokb_dir/allAnnotatedVariants.tsv

wget https://www.oncokb.org/api/v1/utils/allAnnotatedVariants.tsv -O oncokb_variants.tsv

# Cancer Hotspots (pre-processed)
cancer_hotspots_dir="$PWD/data/reference/cancer_hotspots"
mkdir -p ${cancer_hotspots_dir}
cd ${cancer_hotspots_dir}
wget https://www.cancerhotspots.org/files/hotspots_v2.xls

# From Google DeepMind (choose GRCh38 or GRCh37)
cancer_hotspots_dir="$PWD/data/reference/AlphaMissense"
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
gunzip AlphaMissense_hg38.tsv.gz