workflow PREPARE_SAMPLE {
    take:
        input_csv  // Path to the input CSV file

    main:
        // Create unified channel from the CSV with support for both input types
        sample_list = Channel.fromPath(input_csv)
            .ifEmpty { exit(1, "Samplesheet not found: ${input_csv}") }
            .splitCsv(header: true)
        
        input_samples = sample_list
            .map { 
                row ->
                // Extract sample info
                def patient_id = row.patient ? row.patient.trim() : null
                def sample_id = row.sample ? row.sample.trim() : null
                def status = row.status ? row.status.trim() : null

                // Validate required fields
                if (!patient_id) error("Missing or empty 'patient' field in row: ${row}")
                if (!sample_id) error("Missing or empty 'sample' field in row: ${row}")
                if (!status) error("Missing or empty 'status' field in row: ${row}")

                // Determine available input types
                def has_fastq = row.containsKey('fastq_1') && 
                                row.containsKey('fastq_2') && 
                                row.fastq_1 && row.fastq_2

                def has_bam = row.containsKey('bam') && row.bam

                // Exit if no input is provided
                if (!has_fastq && !has_bam) {
                    
                    error("Either FASTQ or BAM files must be provided for sample ${sample_id}")
                }
                
                def status_int = status.toInteger()
                
                def meta = [
                    patient_id: patient_id,
                    sample_id: sample_id,
                    status: status_int
                ]
                
                // Initialize all file variables with null
                def fastq_1 = null
                def fastq_2 = null
                def bam = null
                def bai = null

                if (has_fastq) {
                    fastq_1 = file(row.fastq_1)
                    fastq_2 = file(row.fastq_2)
                }
                
                if (has_bam) {
                    bam = file(row.bam)
                    bai = file(row.bai)
                }

                // Return unified format with all possible inputs
                return [meta, fastq_1, fastq_2, bam, bai]
        
            }

        // Prepare the fastq channel for preprocessing and mapping
        fastq = input_samples.
            map { 
                meta, fastq_1, fastq_2, _bam, _bai ->
                    
                    [meta, fastq_1, fastq_2]
                
            }

        tumour_samples = input_samples
            .filter { 
                meta, _fastq_1, _fastq_2, _bam, _bai ->
                meta.status == 1 
            }

        normal_samples = input_samples
            .filter { 
                meta, _fastq_1, _fastq_2, _bam, _bai ->
                meta.status == 0 
            }

        bam_tumour_paired = tumour_samples
            .map {
                meta, fastq_1, fastq_2, bam, bai -> 
                [meta.patient_id, [meta, fastq_1, fastq_2, bam, bai]] 
            }
            .combine(
                normal_samples
                    .map { 
                        meta, fastq_1, fastq_2, bam, bai -> 
                        [meta.patient_id, [meta, fastq_1, fastq_2, bam, bai]] 
                    },
                by: 0
            )
            .map { 
                patient_id, tumour, normal ->

                def tumour_meta = tumour[0]
                def normal_meta = normal[0]

                def meta = [
                    id: "${tumour_meta.sample_id}_vs_${normal_meta.sample_id}",
                    patient_id: patient_id,
                    tumour_id: tumour_meta.sample_id,
                    normal_id: normal_meta.sample_id,
                    is_paired: true
                ]
                
                return [meta, tumour[3], tumour[4], normal[3], normal[4]]
            }

        bam_tumour_paired
            .count()
            .subscribe { count ->
                
                log.info("Found ${count} paired samples with BAM input.")
            }

        bam_tumour_only
        
        // // Create unpaired samples channel
        // normal_patient_ids = normal_samples
        //     .map { it.patient_id }
        //     .unique()
        //     .collect()
        //     .map { ids -> ids.sort() }

        // unpaired_samples = tumour_samples
        //     .branch {
        //         def normal_ids = normal_patient_ids.val
        //         paired: normal_ids.contains(it.patient_id)
        //         unpaired: true
        //     }
        //     .unpaired
        //     .map { tumour -> 
        //         def meta = [
        //             id: tumour.sample_id,
        //             patient_id: tumour.patient_id,
        //             tumour_id: tumour.sample_id,
        //             normal_id: null,
        //             is_paired: false
        //         ]
                
        //         [
        //             meta, 
        //             tumour.bam, tumour.bai, null, null,
        //             tumour.fastq_1, tumour.fastq_2, null, null,
        //             tumour.input_type, null
        //         ]
        //     }
        
        // // Combine paired and unpaired samples
        // all_samples = paired_samples.mix(unpaired_samples)

        // // Create sample channels for variant calling and other downstream analysis
        // // These channels will include both FASTQ-derived and BAM-derived samples
        
        // // For variant callers that need BAM files
        // variant_calling_paired_samples = paired_samples.map { meta, t_bam, t_bai, n_bam, n_bai, 
        //                                                       t_fq1, t_fq2, n_fq1, n_fq2, 
        //                                                       t_type, n_type ->
        //     [meta, t_bam, t_bai, n_bam, n_bai]
        // }
        
        // variant_calling_tumour_only = unpaired_samples.map { meta, t_bam, t_bai, n_bam, n_bai, 
        //                                                      t_fq1, t_fq2, n_fq1, n_fq2, 
        //                                                      t_type, n_type ->
        //     [meta, t_bam, t_bai]
        // }
        
        // // Log sample counts
        // paired_samples.count().subscribe { count ->
        //     log.info("Found ${count} paired samples.")
        // }
        
        // unpaired_samples.count().subscribe { count ->
        //     log.info("Found ${count} unpaired samples.")
        // }
        
        // fastq_samples.count().subscribe { count ->
        //     log.info("Found ${count} samples with FASTQ input (need alignment).")
        // }
        
        // input_samples_by_type.bam.count().subscribe { count ->
        //     log.info("Found ${count} samples with BAM input (skip alignment).")
        // }

    emit:
        // Raw input samples by type for preprocessing
        input_samples            = null
        // fastq_samples            = fastq_samples
        
        // // General sample groupings
        // all_samples              = all_samples
        // tumour_samples           = tumour_samples
        // normal_samples           = normal_samples
        // paired_samples           = paired_samples
        // unpaired_samples         = unpaired_samples
        
        // // Variant calling ready samples
        // variant_calling_paired   = variant_calling_paired_samples
        // variant_calling_unpaired = variant_calling_tumour_only
}

workflow  {
    
    PREPARE_SAMPLE(params.input)
    
}