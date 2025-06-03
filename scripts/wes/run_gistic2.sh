#!/bin/bash
#SBATCH --job-name=GISTIC2_parallel
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
module load gnuparallel/20211222

# Set paths
gistic2_dir=/share1/GISTIC/2.0.23/bin
reference="/lustre1/g/path_my/250224_DFSP_WES/data/reference/hg38.UCSC.add_miR.160920.refgene.mat"
work_dir="/lustre1/g/path_my/250224_DFSP_WES/data/wes/variant_calling/cnv/gistic2"

# Create a simple function to run GISTIC2 on one file
run_gistic() {
    local input_file=$1
    local basename=$(basename "$input_file" .seg)
    local out_dir="${work_dir}/gistic2_${basename}"
    
    echo "Processing: $input_file"
    mkdir -p "$out_dir"
    
    # Run GISTIC2
    ${gistic2_dir}/gistic2 \
        -seg "${input_file}" \
        -b "${out_dir}" \
        -refgene "${reference}" \
        -qvt 0.1 \
        -ta 0.3 \
        -td 0.3 \
        -conf 0.99 \
        &> "${out_dir}/gistic2.log"
    
    echo "Finished: $input_file"
}

# Export the function for GNU Parallel
export -f run_gistic
export gistic2_dir
export reference
export work_dir

# Count the number of seg files
file_count=$(find "$work_dir" -name "*.seg" | wc -l)
echo "Found $file_count .seg files"

# Calculate number of parallel jobs
if [ $file_count -ge 15 ]; then
    jobs=15
else
    jobs=$file_count
fi

jobs=1
echo "Will run with $jobs parallel jobs"

# Find all ready-to-use segment files and run in parallel
find "$work_dir" -name "*.seg" | parallel \
    --jobs ${jobs} \
    --progress \
    run_gistic

echo "All GISTIC2 analyses completed"