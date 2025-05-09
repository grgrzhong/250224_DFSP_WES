workflow PREPARE_SAMPLE {
    take:
        input_csv  // Path to the input CSV file

    main:
        // Create unified channel from the CSV with support for both input types
        input_samples = Channel.fromPath(input_csv)
            .ifEmpty { exit(1, "Samplesheet not found: ${input_csv}") }
            .splitCsv(header: true)
            .map { row ->
                // Extract sample info
                def patient_id = row.patient ? row.patient.trim() : null
                def sample_id = row.sample ? row.sample.trim() : null
                def status = row.status ? row.status.trim() : null

                // Validate required fields
                if (!patient_id) error("Missing or empty 'patient' field in row: ${row}")
                if (!sample_id) error("Missing or empty 'sample' field in row: ${row}")
                if (!status) error("Missing or empty 'status' field in row: ${row}")

                // Determine available input types
                def has_fastq = row.containsKey('fastq_1') && row.containsKey('fastq_2') && 
                               row.fastq_1 && row.fastq_2
                def has_bam = row.containsKey('bam') && row.bam

                // Exit if no input is provided
                if (!has_fastq && !has_bam) {
                    error("Either FASTQ or BAM files must be provided for sample ${sample_id}")
                }

                // Determine which input type to use based on priority
                // Default priority: BAM > FASTQ (unless overridden by params)
                def use_bam_priority = params.containsKey('input_priority') ? 
                                      params.input_priority == 'bam' : true
                                  
                def input_type
                if (has_bam && (use_bam_priority || !has_fastq)) {
                    input_type = "bam"
                    log.info("Using BAM input for sample ${sample_id}")
                } else if (has_fastq) {
                    input_type = "fastq"
                    log.info("Using FASTQ input for sample ${sample_id}")
                }
                
                def status_int = status.toInteger()
                
                def result = [
                    patient_id: patient_id,
                    sample_id: sample_id,
                    status: status_int,
                    input_type: input_type
                ]
                
                // Add files based on selected input type
                if (input_type == "fastq") {
                    def fastq_1 = file(row.fastq_1)
                    def fastq_2 = file(row.fastq_2)
                    
                    if (!fastq_1.exists()) error("FASTQ R1 file not found: ${fastq_1}")
                    if (!fastq_2.exists()) error("FASTQ R2 file not found: ${fastq_2}")
                    
                    result.fastq_1 = fastq_1
                    result.fastq_2 = fastq_2
                    result.bam = null
                    result.bai = null
                } else {
                    def bam = file(row.bam)
                    def bai = row.containsKey('bai') && row.bai ? file(row.bai) : file("${row.bam}.bai")
                    
                    if (!bam.exists()) error("BAM file not found: ${bam}")
                    if (!bai.exists()) error("BAI file not found: ${bai}")
                    
                    result.fastq_1 = null
                    result.fastq_2 = null
                    result.bam = bam
                    result.bai = bai
                }
                
                return result
            }

         // Split by input type
        input_samples_by_type = input_samples.branch {
            fastq: it.input_type == "fastq"
            bam: it.input_type == "bam"
        }
        
        // Extract FASTQ samples specifically for preprocessing
        fastq_samples = input_samples_by_type.fastq
            .map { sample -> 
                def meta = [
                    id: sample.sample_id, 
                    patient_id: sample.patient_id,
                    status: sample.status
                ]
                [meta, sample.fastq_1, sample.fastq_2]
            }
            
        // Split samples into tumour and normal
        tumour_samples = input_samples.filter { it.status == 1 }
        normal_samples = input_samples.filter { it.status == 0 }

        // Create paired samples channel
        paired_samples = tumour_samples
            .map { tumour -> [tumour.patient_id, tumour] }
            .combine(
                normal_samples.map { normal -> [normal.patient_id, normal] },
                by: 0
            )
            .map { patient_id, tumour, normal -> 
                def meta = [
                    id: "${tumour.sample_id}_vs_${normal.sample_id}",
                    patient_id: patient_id,
                    tumour_id: tumour.sample_id,
                    normal_id: normal.sample_id,
                    is_paired: true
                ]
                
                [
                    meta, 
                    tumour.bam, tumour.bai, normal.bam, normal.bai,
                    tumour.fastq_1, tumour.fastq_2, normal.fastq_1, normal.fastq_2,
                    tumour.input_type, normal.input_type
                ]
            }
        
        // Create unpaired samples channel
        normal_patient_ids = normal_samples
            .map { it.patient_id }
            .unique()
            .collect()
            .map { ids -> ids.sort() }

        unpaired_samples = tumour_samples
            .branch {
                def normal_ids = normal_patient_ids.val
                paired: normal_ids.contains(it.patient_id)
                unpaired: true
            }
            .unpaired
            .map { tumour -> 
                def meta = [
                    id: tumour.sample_id,
                    patient_id: tumour.patient_id,
                    tumour_id: tumour.sample_id,
                    normal_id: null,
                    is_paired: false
                ]
                
                [
                    meta, 
                    tumour.bam, tumour.bai, null, null,
                    tumour.fastq_1, tumour.fastq_2, null, null,
                    tumour.input_type, null
                ]
            }
        
        // Combine paired and unpaired samples
        all_samples = paired_samples.mix(unpaired_samples)

        // Create sample channels for variant calling and other downstream analysis
        // These channels will include both FASTQ-derived and BAM-derived samples
        
        // For variant callers that need BAM files
        variant_calling_paired_samples = paired_samples.map { meta, t_bam, t_bai, n_bam, n_bai, 
                                                              t_fq1, t_fq2, n_fq1, n_fq2, 
                                                              t_type, n_type ->
            [meta, t_bam, t_bai, n_bam, n_bai]
        }
        
        variant_calling_tumour_only = unpaired_samples.map { meta, t_bam, t_bai, n_bam, n_bai, 
                                                             t_fq1, t_fq2, n_fq1, n_fq2, 
                                                             t_type, n_type ->
            [meta, t_bam, t_bai]
        }
        
        // Log sample counts
        paired_samples.count().subscribe { count ->
            log.info("Found ${count} paired samples.")
        }
        
        unpaired_samples.count().subscribe { count ->
            log.info("Found ${count} unpaired samples.")
        }
        
        fastq_samples.count().subscribe { count ->
            log.info("Found ${count} samples with FASTQ input (need alignment).")
        }
        
        input_samples_by_type.bam.count().subscribe { count ->
            log.info("Found ${count} samples with BAM input (skip alignment).")
        }

    emit:
        // Raw input samples by type for preprocessing
        fastq_samples            = fastq_samples
        
        // General sample groupings
        all_samples              = all_samples
        tumour_samples           = tumour_samples
        normal_samples           = normal_samples
        paired_samples           = paired_samples
        unpaired_samples         = unpaired_samples
        
        // Variant calling ready samples
        variant_calling_paired   = variant_calling_paired_samples
        variant_calling_unpaired = variant_calling_tumour_only
}