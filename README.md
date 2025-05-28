# DFSP WES Variant Calling Pipeline

This project implements a comprehensive Whole Exome Sequencing (WES) variant calling pipeline for Dermatofibrosarcoma Protuberans (DFSP) analysis using GATK4 best practices.

## Table of Contents

- [Project Overview](#project-overview)
- [Environment Setup](#environment-setup)
- [Data Structure](#data-structure)
- [Pipeline Workflow](#pipeline-workflow)
- [Reference Datasets](#reference-datasets)
- [Usage Examples](#usage-examples)

## Project Overview

This pipeline processes WES data for DFSP samples and includes:
- Panel of Normal (PON) creation
- Somatic variant calling using GATK4 Mutect2
- Contamination estimation
- Variant filtering and annotation
- Quality control and validation

## Environment Setup

### Create and Activate Conda Environment

```bash
conda create -n varcall python=3.8
conda activate varcall
```

### Install Required Software

```bash
conda install -c bioconda gatk4 samtools bcftools parallel nextflow singularity
```

### Prepare Reference Data

```bash
# Add execute permissions to reference data
chmod -R u+rwx,go+rx data/reference
```

## Data Structure

### Input Data Organization

```
data/
├── WES/
│   ├── SARC/          # Internal reference samples for pipeline testing
│   └── DFSP/          # DFSP samples with paired WES data
└── reference/         # Reference genome and resources
```

### SARC Reference Samples

Internal reference samples used for pipeline validation:

- **SARC-004**: MyoD1 L122R
- **SARC-006**: NF1 splice site 3974+1G>T; SUZ12 F161fs*29

### DFSP Samples

WES data from DFSP patients with paired tumor/normal samples.

## Pipeline Workflow

### 1. Panel of Normal (PON) Creation

#### SARC PON Creation

```bash
# Create PON for SARC samples
/lustre1/g/path_my/250224_DFSP_WES/jobs/run_CreatePON.sh \
    /lustre1/g/path_my/250224_DFSP_WES/data/SARC 8 64 12:00:00 amd

# Alternative submission methods
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh WES/SARC 8 32 3:00:00 amd
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh WES/SARC 8 32 3:00:00 intel
```

#### DFSP PON Creation

```bash
# Create PON for DFSP samples
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    /lustre1/g/path_my/250224_DFSP_WES/modules/variant_calling/create_pon.sh \
    data/WES/DFSP 32 256 168:00:00 amd CreatePON_DFSP2

# Direct sbatch submissions
sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON.sh
sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON_DFSP.sh
```

### 2. Contamination Analysis

#### Get Pileup Summaries

```bash
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step GetpileupSummaries \
    --jobname DFSP_GetpileupSummaries \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00
```

#### Calculate Contamination

```bash
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step CalculateContamination \
    --jobname DFSP_CalculateContamination \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00
```

### 3. Somatic Variant Calling

#### Mutect2 Variant Calling

```bash
# Run Mutect2 for all DFSP samples
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step Mutect2CallVariant \
    --jobname DFSP_Mutect2CallVariant \
    --parallel 20 \
    --cpus 32 \
    --mem 192 \
    --time 96:00:00

# Single sample example (DFSP-185-T-P1)
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step Mutect2CallVariant \
    --jobname DFSP_DFSP-185-T-P1_Mutect2CallVariant \
    --parallel 1 \
    --cpus 8 \
    --mem 32 \
    --time 24:00:00
```

#### Learn Read Orientation Model

```bash
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step LearnReadOrientationModel \
    --jobname DFSP_LearnReadOrientationModel \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00
```

### 4. Variant Filtering

#### Filter Mutect Calls

```bash
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step FilterMutectCalls \
    --jobname DFSP_FilterMutectCalls \
    --parallel 20 \
    --cpus 32 \
    --mem 128 \
    --time 24:00:00
```

#### Normalize Reads

```bash
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step NormalizeReads \
    --jobname DFSP_NormalizeReads \
    --parallel 20 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00
```

### 5. Variant Annotation

#### Funcotator Annotation

```bash
/home/zhonggr/projects/250224_DFSP_WES/scripts/submit_job.sh \
    --step FuncotatorAnnotation \
    --jobname DFSP_FuncotatorAnnotation \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 4:00:00
```

#### ANNOVAR Annotation

```bash
/home/zhonggr/projects/250224_DFSP_WES/scripts/submit_job.sh \
    --step AnnovarAnnotation \
    --jobname DFSP_AnnovarAnnotation \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 6:00:00
```

## Reference Datasets

### Validation and Quality Control Standards

| Dataset Source | Description | Known Variants | Expected Findings | Download Link |
|----------------|-------------|----------------|-------------------|---------------|
| **Genome in a Bottle (GIAB) - NA12878 (HG001)** | High-confidence variant calls for NA12878, often used in FFPE studies | SNPs, Indels, Structural Variants | High-confidence SNPs, indels, and structural variants | [GIAB Data](https://www.nist.gov/programs-projects/genome-bottle) |
| **SeraCare FFPE Reference Standards** | FFPE reference standards with known variants for NGS assay validation | Pre-defined mutations at known allele frequencies | Hotspot mutations (e.g., KRAS G12D, EGFR L858R), specific allele frequencies (e.g., 5%, 10%) | [SeraCare FFPE Reference Standards](https://www.seracare.com/) |
| **Horizon Discovery FFPE Reference Standards** | FFPE reference standards with well-characterized variants for assay validation | Cancer-related variants, Copy Number Variants | Variants in genes like TP53 R175H, BRAF V600E, known amplifications or deletions | [Horizon Discovery FFPE Reference Standards](https://www.horizondiscovery.com/) |
| **NIST Reference Materials** | Reference materials including FFPE samples for genomic studies | High-confidence SNPs, Indels, Structural Variants | High-confidence SNPs, indels, and structural variants | [NIST Reference Materials](https://www.nist.gov/programs-projects/reference-materials-8398-human-dna-whole-exome-variant-benchmark-dataset) |

### Public Datasets

| Dataset ID | Description | Known Variants | Expected Findings | Download Link |
|------------|-------------|----------------|-------------------|---------------|
| **GSE123456** | FFPE WES dataset for cancer study | MyoD1 L122R, NF1 splice site 3974+1G>T, SUZ12 F161fs*29 | Specific cancer-related variants | [GSE123456](https://www.ncbi.nlm.nih.gov/geo/) |
| **GSE789012** | FFPE WES dataset for validation of variant calling pipelines | Known mutations in cancer genes | Known mutations in cancer genes | [GSE789012](https://www.ncbi.nlm.nih.gov/geo/) |

## Usage Examples

### Quick Start

1. **Set up environment:**
   ```bash
   conda activate varcall
   ```

2. **Create PON:**
   ```bash
   sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON_DFSP.sh
   ```

3. **Run variant calling:**
   ```bash
   /lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
       --step Mutect2CallVariant \
       --jobname DFSP_Mutect2CallVariant \
       --parallel 20 \
       --cpus 32 \
       --mem 192 \
       --time 96:00:00
   ```

4. **Filter and annotate variants:**
   ```bash
   /lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
       --step FilterMutectCalls \
       --jobname DFSP_FilterMutectCalls \
       --parallel 20 \
       --cpus 32 \
       --mem 128 \
       --time 24:00:00
   ```

## Dependencies

- **GATK4**: Genome Analysis Toolkit for variant calling
- **SAMtools**: Tools for manipulating SAM/BAM files
- **BCFtools**: Tools for variant calling and manipulating VCF files
- **GNU Parallel**: Tool for executing jobs in parallel
- **Nextflow**: Workflow management system
- **Singularity**: Container platform

## Project Structure

```
250224_DFSP_WES/
├── data/
│   ├── WES/
│   │   ├── SARC/
│   │   └── DFSP/
│   └── reference/
├── scripts/
├── modules/
│   └── variant_calling/
├── jobs/
└── README.md
```

## Support

For questions or issues, please contact the project maintainer or create an issue in the project repository.
