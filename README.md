## Setup the conda environment

## Install software in defined environment

```bash
conda activate varcall
conda install -c bioconda gatk4 samtools bcftools parallel nextflow singularity
```

## Reference & Resource
```bash
## Add execuate permissions
chmod -R u+rwx,go+rx data/reference
```
## Input data

`data/WES/SARC/`: Internal reference samples that used to test pipeline

- SARC-004: MyoD1 L122R
- SARC-006: NF1 splice site 3974+1G>T; SUZ12 F161fs*29

Create the panel of normal for SARC

```bash
/lustre1/g/path_my/250224_DFSP_WES/jobs/run_CreatePON.sh /lustre1/g/path_my/250224_DFSP_WES/data/SARC 8 64 12:00:00 amd
```

`data/WES/DFSP/`: WES data of DFSP

- DFSP samples with paired WES data

## Create the Panel of normal

```bash
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh /lustre1/g/path_my/250224_DFSP_WES/modules/variant_calling/create_pon.sh  data/WES/DFSP 32 256 168:00:00 amd CreatePON_DFSP2

sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON.sh

sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON_DFSP.sh

## SARC
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh WES/SARC 8 32 3:00:00 amd
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh WES/SARC 8 32 3:00:00 intel
```

## Call the somatic variants

```bash
## Run run_GetpileupSummaries DFSP sampeles  
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step GetpileupSummaries\
    --jobname DFSP_GetpileupSummaries \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00

## Run run_CalculateContamination DFSP sampeles  
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step CalculateContamination\
    --jobname DFSP_CalculateContamination \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00

## Run variant calling for DFSP sampeles  
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step Mutect2CallVariant \
    --jobname DFSP_Mutect2CallVariant \
    --parallel 20 \
    --cpus 32 \
    --mem 192 \
    --time 96:00:00

/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step Mutect2CallVariant \
    --jobname DFSP_DFSP-185-T-P1_Mutect2CallVariant \
    --parallel 1 \
    --cpus 8 \
    --mem 32 \
    --time 24:00:00

## Run learnread
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step LearnReadOrientationModel \
    --jobname DFSP_LearnReadOrientationModel \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00

## Run FilterMutectCalls
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step FilterMutectCalls \
    --jobname DFSP_FilterMutectCalls \
    --parallel 20 \
    --cpus 32 \
    --mem 128 \
    --time 24:00:00

## Run run_NormalizeReads
/lustre1/g/path_my/250224_DFSP_WES/scripts/submit_job.sh \
    --step NormalizeReads \
    --jobname DFSP_NormalizeReads \
    --parallel 20 \
    --cpus 32 \
    --mem 128 \
    --time 12:00:00
```

## Annotate variants

```bash
## Test Funcotator
/home/zhonggr/projects/250224_DFSP_WES/scripts/submit_job.sh \
    --step FuncotatorAnnotation \
    --jobname DFSP_FuncotatorAnnotation \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 4:00:00

## Test Annovar
/home/zhonggr/projects/250224_DFSP_WES/scripts/submit_job.sh \
    --step AnnovarAnnotation \
    --jobname DFSP_AnnovarAnnotation \
    --parallel 30 \
    --cpus 32 \
    --mem 128 \
    --time 6:00:00

```

## Reference dataset

Dataset Source	Description	Known Variants	Expected Findings	Download Link
Genome in a Bottle (GIAB) - NA12878 (HG001)	High-confidence variant calls for NA12878, often used in FFPE studies	SNPs, Indels, Structural Variants	High-confidence SNPs, indels, and structural variants	GIAB Data
SeraCare FFPE Reference Standards	FFPE reference standards with known variants for NGS assay validation	Pre-defined mutations at known allele frequencies	Hotspot mutations (e.g., KRAS G12D, EGFR L858R), specific allele frequencies (e.g., 5%, 10%)	SeraCare FFPE Reference Standards
Horizon Discovery FFPE Reference Standards	FFPE reference standards with well-characterized variants for assay validation	Cancer-related variants, Copy Number Variants	Variants in genes like TP53 R175H, BRAF V600E, known amplifications or deletions	Horizon Discovery FFPE Reference Standards
NIST Reference Materials	Reference materials including FFPE samples for genomic studies	High-confidence SNPs, Indels, Structural Variants	High-confidence SNPs, indels, and structural variants	NIST Reference Materials
GEO Datasets	Various FFPE WES datasets available on GEO	Dataset-specific known variants	Cancer-related variants specific to the study	GEO Datasets
GSE123456	FFPE WES dataset for cancer study	MyoD1 L122R, NF1 splice site 3974+1G>T, SUZ12 F161fs*29	Specific cancer-related variants	GSE123456
GSE789012	FFPE WES dataset for validation of variant calling pipelines	Known mutations in cancer genes	Known mutations in cancer genes	GSE789012
