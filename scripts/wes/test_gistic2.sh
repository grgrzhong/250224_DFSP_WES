#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

module load GISTIC/2.0.23

# Set gistic2 path
gistic2_dir=/share1/GISTIC/2.0.23/bin
reference="/lustre1/g/path_my/250224_DFSP_WES/data/reference/hg38.UCSC.add_miR.160920.refgene.mat"
segment_file="/lustre1/g/path_my/250224_DFSP_WES/data/wes/variant_calling/cnv/DFSP_GISTIC_all_tumour.csv"
out_dir="/lustre1/g/path_my/250224_DFSP_WES/data/wes/variant_calling/cnv/gistic2"
mkdir -p "$out_dir"

head -n 5 ${segment_file}

# Fix the file format (run these commands) for GISTIC2
sed 's/"//g' ${segment_file} | tr ',' '\t' > ${segment_file%.*}_fixed.seg

# Correct column order and names (GISTIC2 is case-sensitive!)
# Gistic2 requires the following columns and case-sensitive:
# Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
awk -F '\t' 'BEGIN {OFS="\t"} NR==1 {print "Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean"; next} {print $2, $3, $4, $5, $6, $7}' ${segment_file%.*}_fixed.seg > ${segment_file%.*}_gistic-ready.seg
head -n 5 ${segment_file%.*}_gistic-ready.seg

segment_file="/lustre1/g/path_my/250224_DFSP_WES/data/wes/variant_calling/cnv/DFSP_GISTIC_all_tumour_gistic-ready.seg"

# Run GISTIC2, https://github.com/bzhanglab/GISTIC2_example
# may need to try different parameters
${gistic2_dir}/gistic2 \
    -seg ${segment_file} \
    -b ${out_dir} \
    -refgene ${reference} \
    -qvt 0.1 \
    -ta 0.3 \
    -td 0.3 \
    -conf 0.99 \
    &> ${out_dir}/gistic2.log

# $gistic2_dir/gistic2 \
#     -seg ${segment_file} \
#     -b ${out_dir} \
#     -refgene ${reference} \
#     -genegistic 1 \
#     -smallmem 0 \
#     -rx 0 \
#     -broad 1 \
#     -brlen 0.7 \
#     -conf 0.99 \
#     -armpeel 1 \
#     -savegene 1 \
# 	-gcm extreme \
#     -v 30 \
#     -maxseg 46000 \
#     -ta 0.3 \
#     -td 0.3 \
#     -cap 1.5 \
#     -js 4