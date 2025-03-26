
## General Practise in coding for modules
1. Environment Variables and Constants
Use UPPERCASE for:

 - Environment variables
 - Constants
 - Variables exported to child processes

```bash
# Environment and constant variables
REFERENCE_DIR="/lustre1/g/path_my/Reference"
MAX_THREADS=32
GENOME_VERSION="hg38"
```
2. Local Variables
Use lowercase or snake_case for:

 - Function-local variables
 - Script-local variables
 - Loop variables

```bash
# Local variables
local sample_name="SARC-001"
local tumor_bam="${bam_dir}/${sample_name}.bam"
read_count=0
```
3. Function Names
Use snake_case for:

 - Function names
 - Script names

```bash
annotate_with_funcotator() {
    local sample_id=$1
    local output_dir=$2
    # ...
}
```

4. Example Refactoring of  Code

```bash
#!/bin/bash

# Constants and environment variables
export REFERENCE_DIR="/lustre1/g/path_my/Reference"
export BASE_DIR="/lustre1/g/path_my"
export PROJECT_DIR="250224_DFSP_WES"

# Derived paths (still uppercase as they're exported)
export REFERENCE="${REFERENCE_DIR}/Gencode/gencode.hg38.v36.primary_assembly.fa"
export OUTPUT_DIR="${BASE_DIR}/${PROJECT_DIR}/${WORK_DIR}/Mutect-Call"
export LOG_DIR="${OUTPUT_DIR}/logs"

# Local variables
local avail_mem
if [ -n "$SLURM_MEM_PER_NODE" ]; then
    avail_mem=$((SLURM_MEM_PER_NODE * 80 / 102400))
else
    sys_mem=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    avail_mem=$((sys_mem * 80 / 1024 / 1024))
fi

# Function using local variables
annotate_with_funcotator() {
    local case_id=$1
    local tumor_sample=$2
    local output_dir="${OUTPUT_DIR}/${tumor_sample}"
    
    log_message "Annotating variants for sample: ${tumor_sample}"
    # ...
}
```
5. Summary

UPPERCASE:
 - Environment variables (PATH, HOME)
 - Constants (MAX_THREADS, REFERENCE_DIR)
 - Exported variables (export BAM_DIR)

lowercase/snake_case:
 - Local variables (sample_name, tumor_bam)
 - Function parameters (case_id, output_dir)
 - Loop variables (i, file_name)

Descriptive Names:
 - Use clear, descriptive names
 - Avoid abbreviations unless common
 - Include units if applicable (avail_mem_gb)

Scope Indication:
 - Prefix temporary variables with _
 - Use local keyword in functions
 - Document exported variables

Consistency:
 - Maintain consistent style throughout project
 - Follow team/project conventions
 - Document naming conventions in README

## Optimal Memory Setup for GATK Mutect2 in HPC Environment
When running GATK Mutect2 on an HPC system with GNU Parallel for 20 samples simultaneously, memory allocation requires careful consideration to balance performance and resource usage.

GATK Mutect2 Memory Requirements
- Minimum: 4-6GB per Mutect2 instance
- Recommended: 8-16GB per instance for standard WES samples
- For WGS: 16-32GB per instance may be needed

Factors Affecting Memory Usage
- Sample complexity: Tumor samples with high heterogeneity require more memory
- Depth of coverage: Higher coverage samples need more memory
- Region size: Processing larger genomic regions requires more memory
- Reference genome: hg38 requires more memory than hg19

Optimal `-Xmx` Settings for HPC Environment
Formula for Memory Allocation
For 20 parallel samples, a good approach is:
```bash
avail_mem=$(( (total_node_memory * 0.85) / 20 ))
```

