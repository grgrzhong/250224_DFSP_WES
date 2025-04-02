#!/usr/bin/env nextflow

// Default parameters
params.samplesheet = "${launchDir}/samplesheet.csv"
params.outdir = "${launchDir}/results/sv_calling"
params.reference = "/path/to/reference/genome.fa"
params.exclude = "/path/to/delly/excludeTemplates/human.hg38.excl.tsv"
params.regions = "/path/to/regions.bed.gz"
params.publishMode = "copy"

// Process 1: Call SVs with Delly
process DELLY_CALL {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}/Delly", mode: params.publishMode
    
    input:
    tuple val(patient_id), val(tumor_id), val(normal_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(reference)
    path(exclude)
    
    output:
    tuple val(patient_id), val(tumor_id), val(normal_id), path("${tumor_id}.bcf"), path("samples.tsv"), emit: delly_out
    
    script:
    """
    # Call SVs using Delly
    delly call \\
    -x ${exclude} \\
    -o ${tumor_id}.bcf \\
    -g ${reference} \\
    ${tumor_bam} \\
    ${normal_bam}
    
    # Create samples file for filtering
    echo -e "${tumor_id}\\ttumor\\n${normal_id}\\tcontrol" > samples.tsv
    """
}

// Process 2: Filter Delly results
process DELLY_FILTER {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}/Delly", mode: params.publishMode
    
    input:
    tuple val(patient_id), val(tumor_id), val(normal_id), path(delly_bcf), path(samples_tsv)
    
    output:
    tuple val(patient_id), val(tumor_id), val(normal_id), path("${tumor_id}.somaticSV.delly.vcf.gz"), path("${tumor_id}.somaticSV.delly.vcf.gz.tbi"), emit: delly_vcf
    
    script:
    """
    # Filter for somatic SVs
    delly filter \\
    -f somatic \\
    -o ${tumor_id}.somatic.filter.bcf \\
    -s ${samples_tsv} \\
    ${delly_bcf}
    
    # Convert to VCF and compress
    bcftools view ${tumor_id}.somatic.filter.bcf -Oz > ${tumor_id}.somaticSV.delly.vcf.gz
    
    # Index the VCF
    tabix ${tumor_id}.somaticSV.delly.vcf.gz
    """
}

// Process 3: Call SVs with Manta
process MANTA_CALL {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}/Manta", mode: params.publishMode
    
    input:
    tuple val(patient_id), val(tumor_id), val(normal_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(reference)
    path(regions)
    
    output:
    tuple val(patient_id), val(tumor_id), val(normal_id), path("runWorkflow.py"), path("results"), emit: manta_results
    
    script:
    """
    # Configure Manta
    configManta.py \\
    --normalBam=${normal_bam} \\
    --tumorBam=${tumor_bam} \\
    --exome \\
    --referenceFasta=${reference} \\
    --runDir=. \\
    --callRegions=${regions}
    
    # Run Manta workflow
    python runWorkflow.py -j ${task.cpus}
    """
}

// Process 4: Process and rename Manta outputs
process MANTA_PROCESS {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}/Manta", mode: params.publishMode
    
    input:
    tuple val(patient_id), val(tumor_id), val(normal_id), path(workflow), path(results)
    path(reference)
    
    output:
    tuple val(patient_id), val(tumor_id), val(normal_id), path("${tumor_id}.somaticSV.manta.vcf.gz"), path("${tumor_id}.somaticSV.manta.vcf.gz.tbi"), emit: manta_vcf
    path("${tumor_id}.*.vcf.gz"), emit: manta_other_vcfs
    path("${tumor_id}.*.vcf.gz.tbi"), emit: manta_other_indices
    
    script:
    """
    # Convert inversions
    convertInversion.py \$(which samtools) ${reference} results/variants/somaticSV.vcf.gz
    tabix -f --preset vcf results/variants/somaticSV.vcf.gz
    
    # Rename output files
    cp results/variants/candidateSmallIndels.vcf.gz ${tumor_id}.candidateSmallIndels.vcf.gz
    cp results/variants/candidateSmallIndels.vcf.gz.tbi ${tumor_id}.candidateSmallIndels.vcf.gz.tbi
    cp results/variants/candidateSV.vcf.gz ${tumor_id}.candidateSV.vcf.gz
    cp results/variants/candidateSV.vcf.gz.tbi ${tumor_id}.candidateSV.vcf.gz.tbi
    cp results/variants/diploidSV.vcf.gz ${tumor_id}.diploidSV.vcf.gz
    cp results/variants/diploidSV.vcf.gz.tbi ${tumor_id}.diploidSV.vcf.gz.tbi
    cp results/variants/somaticSV.vcf.gz ${tumor_id}.somaticSV.manta.vcf.gz
    cp results/variants/somaticSV.vcf.gz.tbi ${tumor_id}.somaticSV.manta.vcf.gz.tbi
    
    # Prepare Manta output for merging
    bcftools view \\
    --samples ${tumor_id},${normal_id} \\
    --output-type z \\
    --output-file ${tumor_id}.manta.swap.vcf.gz \\
    ${tumor_id}.somaticSV.manta.vcf.gz
    
    tabix --preset vcf ${tumor_id}.manta.swap.vcf.gz
    """
}

// Process 5: Merge and filter Delly and Manta results
process MERGE_FILTER_SV {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}", mode: params.publishMode
    
    input:
    tuple val(patient_id), val(tumor_id), val(normal_id), path(delly_vcf), path(delly_idx), path(manta_vcf), path(manta_idx)
    
    output:
    tuple val(patient_id), val(tumor_id), val(normal_id), path("${tumor_id}.delly.manta.vcf.gz"), path("${tumor_id}.delly.manta.vcf.gz.tbi"), emit: merged_vcf
    
    script:
    """
    # Concatenate Delly and Manta results
    bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${tumor_id}.delly.manta.unfiltered.vcf.gz \\
    ${delly_vcf} \\
    ${manta_vcf}
    
    tabix --preset vcf ${tumor_id}.delly.manta.unfiltered.vcf.gz
    
    # Filter for PASS variants and sort
    bcftools filter \\
    --include 'FILTER="PASS"' \\
    ${tumor_id}.delly.manta.unfiltered.vcf.gz | \\
    bcftools sort \\
    --output-type z \\
    --output-file ${tumor_id}.delly.manta.vcf.gz
    
    tabix --preset vcf ${tumor_id}.delly.manta.vcf.gz
    """
}

// Process 6: Annotate SVs with AnnotSV
process ANNOTSV {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}", mode: params.publishMode
    
    input:
    tuple val(patient_id), val(tumor_id), val(normal_id), path(merged_vcf), path(merged_idx)
    
    output:
    path("${tumor_id}.somaticSV.annotated.tsv"), emit: annotated_sv
    path("${tumor_id}.AnnotSV.log"), emit: log
    
    script:
    """
    # Annotate SVs
    AnnotSV \\
    -SVinputFile ${merged_vcf} \\
    -outputFile ${tumor_id}.somaticSV.annotated.tsv \\
    -genomeBuild GRCh38 \\
    -annotationMode both \\
    -SVminSize 50 \\
    >& ${tumor_id}.AnnotSV.log
    """
}

// Define workflow
workflow {
    // Get reference files
    reference = file(params.reference)
    exclude = file(params.exclude)
    regions = file(params.regions)
    
    // Read samplesheet and prepare tumor-normal pairs
    Channel
        .fromPath(params.samplesheet)
        .ifEmpty { exit 1, "Samplesheet not found: ${params.samplesheet}" }
        .splitCsv(header: true)
        .map { row -> 
            def patient_id = row.patient
            def sample_id = row.sample
            def status = row.status
            def bam = file(row.bam)
            def bai = file(row.bai)
            
            return [patient_id, sample_id, status, bam, bai]
        }
        .branch {
            tumor: it[2] == '1'
            normal: it[2] == '0'
        }
        .set { samples_ch }

    // Match tumor and normal samples by patient ID to create tumor-normal pairs
    tumor_normal_pairs = samples_ch.tumor
        .cross(samples_ch.normal) { it[0] } // Cross by patient_id
        .map { tumor, normal ->
            def patient_id = tumor[0]
            def tumor_id = tumor[1]
            def tumor_bam = tumor[3]
            def tumor_bai = tumor[4]
            def normal_id = normal[1]
            def normal_bam = normal[3]
            def normal_bai = normal[4]
            
            return [patient_id, tumor_id, normal_id, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    // Step 1: Call SVs with Delly
    DELLY_CALL(tumor_normal_pairs, reference, exclude)
    
    // Step 2: Filter Delly results
    DELLY_FILTER(DELLY_CALL.out.delly_out)
    
    // Step 3: Call SVs with Manta
    MANTA_CALL(tumor_normal_pairs, reference, regions)
    
    // Step 4: Process and rename Manta outputs
    MANTA_PROCESS(MANTA_CALL.out.manta_results, reference)
    
    // Step 5: Combine Delly and Manta outputs
    delly_manta_inputs = DELLY_FILTER.out.delly_vcf
        .join(MANTA_PROCESS.out.manta_vcf, by: [0, 1, 2]) // Join by patient_id, tumor_id, normal_id
    
    // Step 6: Merge and filter
    MERGE_FILTER_SV(delly_manta_inputs)
    
    // Step 7: Annotate SVs
    ANNOTSV(MERGE_FILTER_SV.out.merged_vcf)
}