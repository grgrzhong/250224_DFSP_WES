#!/bin/bash

####
# Download the VEP cache file

vep_data_dir="$PWD/data/Reference/ensembl_vep"
mkdir -p vep_data_dir
cd ${vep_data_dir}
wget -O http://ftp.ensembl.org/pub/release-113/variation/vep/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf *.tar.gz
rm homo_sapiens_vep_110_GRCh38.tar.gz # Optionally, remove the tar.gz file after extraction
echo "VEP cache downloaded and extracted to $VEP_DATA_DIR"




