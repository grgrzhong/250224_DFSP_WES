#!/bin/bash

#############################################################################
# Pipeline Step Submission Script
# Description: Run one or more steps of the variant discovery pipeline in HPC
# authors: Zhong Guorui
# date: 2024-03-12
#############################################################################

# Load the project configuration first
PROJECT_CONFIG="/lustre1/g/path_my/250224_DFSP_WES/conf/project_config.sh"
if [ ! -f "$PROJECT_CONFIG" ]; then
    echo "Error: Project configuration not found at $PROJECT_CONFIG"
    exit 1
fi

source "$PROJECT_CONFIG"

# Fixed paths (based on project configuration)
readonly LIB_SCRIPT="${BASE_DIR}/${PROJECT_DIR}/scripts/mutation_calling.sh"
readonly CONDA_INIT="/home/zhonggr/miniforge3/etc/profile.d/conda.sh"
readonly CONDA_ENV="varcall"
readonly SLURM_LOG_DIR="${BASE_DIR}/${PROJECT_DIR}/logs"

# Default SLURM settings - can be overridden via command line
STEP=""
NODES=1
NTASKS=1
CPUS=32
MEM=128
TIME="24:00:00"
PARTITION="amd"
EMAIL="zhonggr@hku.hk"
DEPENDENCY=""
JOB_NAME=""  # Custom job name
PARALLEL=""  # Custom parallel jobs value

# Create log directory
mkdir -p "$SLURM_LOG_DIR"

# Display help
show_help() {
    # Get available steps from mutation_calling.sh
    AVAILABLE_STEPS=$(grep -o "^run_[A-Za-z0-9]*" "$LIB_SCRIPT" | sed 's/^run_//' | sort | tr '\n' ',' | sed 's/,$//')
    
    cat << EOF
Usage: $(basename $0) [options]

Options:
    --step STEP_LIST    Step(s) to run (required)
                        Format: Step1,Step2,Step3 (no spaces)
                        Available: ${AVAILABLE_STEPS}
    --jobname NAME     Custom job name (default: step names)
    --parallel N        Override PARALLEL_JOBS value for GNU parallel
    --nodes N           Number of nodes (default: ${NODES})
    --ntasks N          Number of tasks (default: ${NTASKS})
    --cpus N            Number of CPUs (default: ${CPUS})
    --mem N             Memory in GB (default: ${MEM})
    --time HH:MM:SS     Time limit (default: ${TIME})
    --partition NAME    SLURM partition (default: ${PARTITION})
    --email EMAIL       Notification email (default: ${EMAIL})
    --dependency JOBID  Job dependency (wait for job to complete)
    --help              Show this help message

Examples:
    # Run a single step:
    $(basename $0) --step Mutect2CallVariant

    # Run multiple steps sequentially:
    $(basename $0) --step GetpileupSummaries,CalculateContamination,Mutect2CallVariant

    # Run with custom job name:
    $(basename $0) --step Mutect2CallVariant --jobname SARC_Mutect2
    
    # Run with custom resources:
    $(basename $0) --step GetpileupSummaries --cpus 16 --mem 64
    
    # Run with custom parallel jobs:
    $(basename $0) --step Mutect2CallVariant --parallel 16
    
    # Run after another job completes:
    $(basename $0) --step FilterMutectCalls --dependency 123456
EOF
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --step)       STEP="$2";      shift 2 ;;
            --jobname)   JOB_NAME="$2";   shift 2 ;;
            --parallel)   PARALLEL="$2";   shift 2 ;;
            --cpus)       CPUS="$2";       shift 2 ;;
            --mem)        MEM="$2";        shift 2 ;;
            --time)       TIME="$2";       shift 2 ;;
            --partition)  PARTITION="$2";  shift 2 ;;
            --nodes)      NODES="$2";      shift 2 ;;
            --ntasks)     NTASKS="$2";     shift 2 ;;
            --email)      EMAIL="$2";      shift 2 ;;
            --dependency) DEPENDENCY="$2"; shift 2 ;;
            --help)       show_help;       exit 0  ;;
            *)
                echo "Unknown option: $1"
                show_help
                exit 1
                ;;
        esac
    done

    # Validate required arguments
    if [[ -z "$STEP" ]]; then
        echo "Error: --step is required"
        show_help
        exit 1
    fi
    
    # Verify each step exists in workflow script
    IFS=',' read -ra STEP_ARRAY <<< "$STEP"
    for step in "${STEP_ARRAY[@]}"; do
        if ! grep -q "^run_${step}()" "$LIB_SCRIPT"; then
            echo "Error: Step '${step}' not found in ${LIB_SCRIPT}"
            show_help
            exit 1
        fi
    done
    
    # Validate that parallel is a number if specified
    if [[ -n "$PARALLEL" ]] && ! [[ "$PARALLEL" =~ ^[0-9]+$ ]]; then
        echo "Error: --parallel must be a number"
        exit 1
    fi
}

# Create job script template for multiple steps
get_script_template() {
    local steps=$1
    
    cat << EOF
#!/bin/bash

#############################################################################
# WES Pipeline: Sequential Steps
# Generated: $(date)
#############################################################################

# Load conda environment
if [ -f "${CONDA_INIT}" ]; then
    source "${CONDA_INIT}"
    conda activate ${CONDA_ENV} || { echo "Error: Failed to activate conda environment"; exit 1; }
else
    echo "Error: Conda initialization file not found"
    exit 1
fi

# Source project configuration - This sets all required paths and variables
PROJECT_CONFIG="${PROJECT_CONFIG}"
source "\${PROJECT_CONFIG}"

# Override PARALLEL_JOBS if specified
if [[ -n "${PARALLEL}" ]]; then
    export PARALLEL_JOBS=${PARALLEL}
    echo "Overriding PARALLEL_JOBS value with: ${PARALLEL}"
fi

# # Setup memory usage for GATK java options
# if [[ -n "${MEM}" ]]; then
#     # Calculate available memory per job
#     export AVAIL_MEM=$(($MEM * 80 / 100 / $PARALLEL_JOBS))
#     echo "GATK java memory per job: ${AVAIL_MEM}GB"
    
#     # Calculate maximum memory (80% of total)
#     export MAX_MEM=$(($MEM * 80 / 100))
#     echo "Setting MAX_MEM to: ${MAX_MEM}GB"
# fi

# if [ "$AVAIL_MEM" -lt 4 ]; then
#     export AVAIL_MEM=4 # Set minimum memory threshold
# fi

# Verify TUMOUR_SAMPLES was set from project_config.sh
if [[ -z "\$TUMOUR_SAMPLES" ]]; then
    echo "ERROR: TUMOUR_SAMPLES not defined in project_config.sh"
    exit 1
fi

# Display project configuration
echo "Using project configuration from \${PROJECT_CONFIG}"
echo "Working directory: \${WORK_DIR}"
echo "Tumor samples: \$(echo \$TUMOUR_SAMPLES | wc -w)"
echo "Parallel jobs: \${PARALLEL_JOBS}"

# Source the variant discovery workflow (without executing the main function)
source <(cat "${LIB_SCRIPT}" | sed 's/main "\$@"//g')

# Initialize step tracking variables
declare -A STEP_STATUS
declare -A STEP_RUNTIME

# Run the specified steps sequentially
EOF

    # Add each step to execute sequentially
    IFS=',' read -ra STEP_ARRAY <<< "$steps"
    for step in "${STEP_ARRAY[@]}"; do
        cat << EOF

#############################################################################
# Step: ${step}
#############################################################################
echo "Starting ${step} at \$(date)"
STEP_START_TIME=\$(date +%s)

run_${step}
STEP_STATUS["${step}"]=\$?

STEP_END_TIME=\$(date +%s)
STEP_RUNTIME["${step}"]=\$((\$STEP_END_TIME - \$STEP_START_TIME))

if [ \${STEP_STATUS["${step}"]} -ne 0 ]; then
    echo "ERROR: Step ${step} failed with status \${STEP_STATUS["${step}"]}"
    exit \${STEP_STATUS["${step}"]}
fi
echo "${step} completed successfully at \$(date)"
echo "Runtime: \$(printf '%02dh:%02dm:%02ds' \$((\${STEP_RUNTIME["${step}"]} / 3600)) \$((\${STEP_RUNTIME["${step}"]} % 3600 / 60)) \$((\${STEP_RUNTIME["${step}"]} % 60)))"
echo ""
EOF
    done
    
    # Add step summary report
    cat << EOF

#############################################################################
# Step Execution Summary
#############################################################################
echo "=========== Step Execution Summary ============"
echo "Job start time: $(date)"
echo "Job completion time: \$(date)"
echo ""
echo "Step Status and Runtime:"
echo "-------------------------"
EOF

    # Add each step to the summary report
    for step in "${STEP_ARRAY[@]}"; do
        cat << EOF
echo "${step}: SUCCESS (Runtime: \$(printf '%02dh:%02dm:%02ds' \$((\${STEP_RUNTIME["${step}"]} / 3600)) \$((\${STEP_RUNTIME["${step}"]} % 3600 / 60)) \$((\${STEP_RUNTIME["${step}"]} % 60))))"
EOF
    done
    
    # Add final success message and conda deactivation
    cat << EOF
echo "================================================"
echo "All requested steps completed successfully."

# Deactivate conda environment
conda deactivate

exit 0
EOF
}

# Create job script for steps
create_job_script() {
    local steps=$1
    local step_names=$(echo "$steps" | tr ',' '-')
    
    # Use custom job name if provided, otherwise use step names
    local name_prefix=${JOB_NAME:-${step_names}}
    local job_script="${SLURM_LOG_DIR}/$(date +%Y%m%d)_${name_prefix}_command.sh"
    
    # Create script from template
    get_script_template "$steps" > "${job_script}"
    chmod +x "${job_script}"
    echo "${job_script}"
}

# Submit job to SLURM
submit_job() {
    local job_script=$1
    local steps=$2
    local step_names=$(echo "$steps" | tr ',' '-')
    
    # Use custom job name if provided, otherwise use step names
    local name_prefix=${JOB_NAME:-${step_names}}
    
    # Common SLURM arguments
    local slurm_args=(
        "--job-name=${name_prefix}"
        "--partition=${PARTITION}"
        "--qos=normal"
        "--time=${TIME}"
        "--nodes=${NODES}"
        "--ntasks=${NTASKS}"
        "--cpus-per-task=${CPUS}"
        "--mem=${MEM}G"
        "--mail-type=BEGIN,END,FAIL"
        "--mail-user=${EMAIL}"
        "--output=${SLURM_LOG_DIR}/$(date +%Y%m%d)_${name_prefix}_%j.out"
        "--error=${SLURM_LOG_DIR}/$(date +%Y%m%d)_${name_prefix}_%j.err"
    )
    
    # Add dependency if provided
    if [[ -n "${DEPENDENCY}" ]]; then
        slurm_args+=("--dependency=afterok:${DEPENDENCY}")
    fi
    
    # Submit job and capture job ID
    local job_id=$(sbatch "${slurm_args[@]}" "${job_script}" | awk '{print $4}')
    echo "${job_id}"
}

# Log configuration and execution details
log_execution() {
    local step_count=$(echo "$STEP" | tr ',' '\n' | wc -l)
    
    echo "====================== Job Submission ======================"
    if [[ -n "${JOB_NAME}" ]]; then
        echo "Job name:          ${JOB_NAME}"
    fi
    
    if [[ $step_count -eq 1 ]]; then
        echo "Step:              ${STEP}"
    else
        echo "Steps (${step_count}):"
        IFS=',' read -ra STEP_ARRAY <<< "$STEP"
        for i in "${!STEP_ARRAY[@]}"; do
            echo "  $((i+1)). ${STEP_ARRAY[$i]}"
        done
    fi
    echo "Working directory: ${WORK_DIR}"
    echo "Resources:"
    echo "  CPUs:              ${CPUS}"
    echo "  Memory:            ${MEM} GB"
    echo "  Time:              ${TIME}"
    echo "  Partition:         ${PARTITION}"
    if [[ -n "${PARALLEL}" ]]; then
        echo "  Parallel jobs:     ${PARALLEL} (overridden)"
    else
        echo "  Parallel jobs:     ${PARALLEL_JOBS} (from config)"
    fi
    if [[ -n "${DEPENDENCY}" ]]; then
        echo "  Depends on:      Job ${DEPENDENCY}"
    fi
    echo "Workflow script:     ${LIB_SCRIPT}"
    echo "Log directory:       ${SLURM_LOG_DIR}"
}

# Main function
main() {
    # Check if workflow script exists
    if [[ ! -f "$LIB_SCRIPT" ]]; then
        echo "Error: Library script not found at $LIB_SCRIPT"
        exit 1
    fi
    
    # Parse command line arguments
    parse_args "$@"
    
    # Log execution details
    log_execution
    
    # Create job script
    JOB_SCRIPT=$(create_job_script "${STEP}")
    echo "Job script created: ${JOB_SCRIPT}"
    
    # Submit job
    JOB_ID=$(submit_job "${JOB_SCRIPT}" "${STEP}")
    echo "Job submitted with ID: ${JOB_ID}"
    echo "Check job status with: squeue -j ${JOB_ID}"
    echo "Check the logs at:     ${SLURM_LOG_DIR}"
    echo "==========================================================="
}

# Run main function
main "$@"