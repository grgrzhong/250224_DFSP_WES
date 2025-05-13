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
url="https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1237006675_NR5drEgvZ85edZVWA5an0D0VQBC4&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_rmsk&hgta_ctDesc=table+browser+query+on+rmsk&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbDownBases=200&hgta_doGetBed=get+BED"
wget  -O- -q  "$url" | 
grep -v "#" | 
sort -k1,1 -k2,2n -k3,3n -t$'\t' | 
bgzip -c > repeatmasker.GRCh38.bed.gz
tabix -p bed repeatmasker.GRCh38.bed.gz

## Sort the
gunzip ${ref_dir}/RepeatMasker.bed.gz
# grep -v "#" ${ref_dir}/RepeatMasker.bed

sort -k1,1 -k2,2n -k3,3n -t$'\t' ${ref_dir}/RepeatMasker.bed -o ${ref_dir}/RepeatMasker.sorted.bed
bgzip ${ref_dir}/RepeatMasker.sorted.bed
tabix -p bed ${ref_dir}/RepeatMasker.bed.gz
