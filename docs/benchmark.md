# Benchmark Datasets and Validation Steps

## Benchmark Datasets

| **Category**       | **Dataset**                                | **Description**                                                                 | **Download Link**                                                                 |
|--------------------|--------------------------------------------|---------------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| **SNV/Indel**      | **Genome in a Bottle (GIAB)**              | High-confidence variant calls for benchmarking variant detection pipelines.     | [GIAB](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/)                     |
|                    | **Mills and 1000G Gold Standard Indels**   | A high-confidence set of indels for benchmarking.                              | [Mills Indels](https://ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz) |
|                    | **dbSNP**                                  | A database of known single nucleotide polymorphisms (SNPs) and small indels.    | [dbSNP](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/)                        |
| **Copy Number**    | **Exome Target Regions (BED)**             | BED file defining the regions targeted by your exome capture kit.               | Provided by your capture kit vendor (e.g., Agilent, Illumina, Twist).            |
|                    | **gCNV Population Resources**              | Population-based germline CNV resources for normalization and validation.       | [Broad gCNV Resources](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38) |
| **Structural Variants** | **1000 Genomes SVs**                  | Structural variants identified in the 1000 Genomes Project.                     | [1000 Genomes SVs](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/) |
|                    | **GIAB SVs**                               | High-confidence structural variant calls from Genome in a Bottle.              | [GIAB SVs](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/)                    |
|                    | **DGV (Database of Genomic Variants)**     | A curated database of structural variants from healthy individuals.             | [DGV](http://dgv.tcag.ca/dgv/app/home)                                           |

---

## Validation Steps

### 1. Validation Using Benchmark Datasets

#### a. **SNV/Indel Validation**
- Use high-confidence variant datasets like **Genome in a Bottle (GIAB)** or **Mills and 1000G Gold Standard Indels**.
- **Steps**:
  1. Run your pipeline on the same sample (e.g., NA12878) for which the benchmark dataset is available.
  2. Compare your pipeline's output (VCF file) with the benchmark dataset using tools like:
     - **hap.py**: A tool for comparing VCF files and calculating precision, recall, and F1 score.
       - [hap.py GitHub](https://github.com/Illumina/hap.py)
     - **bcftools isec**: For intersecting and comparing VCF files.
  3. Evaluate metrics like:
     - **Precision**: Fraction of true positives among all called variants.
     - **Recall**: Fraction of true positives among all benchmark variants.
     - **F1 Score**: Harmonic mean of precision and recall.

#### b. **Structural Variant (SV) Validation**
- Use datasets like **GIAB SVs** or **1000 Genomes SVs**.
- **Steps**:
  1. Run your pipeline on the same sample for which the benchmark SV dataset is available.
  2. Compare your SV calls with the benchmark dataset using tools like:
     - **truvari**: A tool for comparing structural variant calls.
       - [truvari GitHub](https://github.com/spiralgenetics/truvari)
  3. Evaluate metrics like:
     - **Precision** and **Recall** for SV detection.
     - Overlap of breakpoints and variant types.

#### c. **Copy Number Variant (CNV) Validation**
- Use datasets like **TCGA CNV Data** or **gCNV Population Resources**.
- **Steps**:
  1. Run your pipeline on the same sample or dataset for which the benchmark CNV dataset is available.
  2. Compare your CNV calls with the benchmark dataset using tools like:
     - **CNVkit**: For CNV analysis and comparison.
       - [CNVkit GitHub](https://github.com/etal/cnvkit)
     - **bedtools intersect**: For comparing CNV regions (BED files).
  3. Evaluate metrics like:
     - Concordance of CNV regions.
     - Precision and recall for CNV detection.

---

### 2. Validation Workflow

#### a. **Prepare Benchmark Datasets**
- Download benchmark datasets for SNVs/indels, SVs, and CNVs.
- Ensure the datasets match your reference genome (e.g., GRCh38 or hg19).

#### b. **Run Your Pipeline**
- Process the same sample(s) as in the benchmark dataset through your pipeline.
- Generate output files (e.g., VCF for SNVs/indels, BED for CNVs, etc.).

#### c. **Compare Results**
- Use comparison tools (e.g., hap.py, truvari, CNVkit) to compare your pipeline's output with the benchmark dataset.
- Generate metrics like precision, recall, and F1 score.

#### d. **Iterate and Optimize**
- If your pipeline's performance is suboptimal, identify the steps causing errors (e.g., alignment, variant calling).
- Optimize parameters or tools in your pipeline and revalidate.

---

### 3. Tools for Validation

| **Tool**       | **Purpose**                          | **Link**                                                                 |
|----------------|--------------------------------------|-------------------------------------------------------------------------|
| **hap.py**     | Compare VCF files for SNV/indel calls | [hap.py GitHub](https://github.com/Illumina/hap.py)                    |
| **truvari**    | Compare structural variant calls     | [truvari GitHub](https://github.com/spiralgenetics/truvari)            |
| **CNVkit**     | CNV analysis and comparison          | [CNVkit GitHub](https://github.com/etal/cnvkit)                        |
| **bedtools**   | Compare genomic regions (e.g., CNVs) | [bedtools GitHub](https://github.com/arq5x/bedtools2)                  |
| **bcftools**   | VCF file manipulation and comparison | [bcftools GitHub](https://github.com/samtools/bcftools)                |

---

### 4. Example Validation Workflow for SNVs/Indels

#### a. **Run Your Pipeline**
```bash
nextflow run  \
    --input /path/to/NA12878_fastq.csv \
    --outdir /path/to/output \
    -profile standard