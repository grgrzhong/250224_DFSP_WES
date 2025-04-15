
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
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh data/WES/DFSP 32 180 120:00:00 amd CreatePON_DFSP

sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON.sh

sbatch /lustre1/g/path_my/250224_DFSP_WES/scripts/CreatePON_DFSP.sh

## SARC
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh WES/SARC 8 32 3:00:00 amd
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh CreatePON.sh WES/SARC 8 32 3:00:00 intel
```

## Call the somatic variants

```bash
## Run variant calling for SARC sampeles  
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh Mutect2_Variant_calling.sh data/WES/SARC 8 64 24:00:00 amd
```

## Annotate variants
```bash
## Annotate variants with GATK Funcotator
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh /lustre1/g/path_my/250224_DFSP_WES/modules/variant_annotation_funcotator.sh data/WES/SARC 8 64 1:00:00 amd 2

## Annotate variants with Annovar
/lustre1/g/path_my/250224_DFSP_WES/submit_job.sh /lustre1/g/path_my/250224_DFSP_WES/modules/variant_annotation_annovar.sh data/WES/SARC 8 64 1:00:00 amd 2
```
