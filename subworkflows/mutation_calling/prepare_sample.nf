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

        input_samples
            .count()
            .subscribe { count ->
                
                log.info("Total samples =  ${count} ")
            }

        // Prepare the fastq channel for preprocessing and mapping
        fastq = input_samples
            .map { 
                meta, fastq_1, fastq_2, _bam, _bai ->
                [meta, fastq_1, fastq_2]
                
            }

        tumour_samples = input_samples
            .filter { 
                meta, _fastq_1, _fastq_2, _bam, _bai ->
                meta.status == 1 
            }

        tumour_samples
            .count()
            .subscribe { count ->
                
                log.info("Total tumour samples = ${count}")
            }

        normal_samples = input_samples
            .filter { 
                meta, _fastq_1, _fastq_2, _bam, _bai ->
                meta.status == 0 
            }
        
        normal_samples
            .count()
            .subscribe { count ->
                
                log.info("Total normal samples = ${count}")
            }

        // Prepare bam input for tumour_normal paired samples
        bam_tumour_normal = tumour_samples
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

        bam_tumour_normal
            .count()
            .subscribe { count ->
                
                log.info("Total tumour-normal samples =  ${count}")
            }

        // Prepare bam input for tumour only samples
        normal_patient_ids = normal_samples
            .map { 
                meta, _fastq_1, _fastq_2, _bam, _bai ->
                meta.patient_id 
            }
            .unique()
            .collect()
            .map { ids -> ids.sort() }

        bam_tumour_only = tumour_samples
            .branch {
                meta, _fastq_1, _fastq_2, _bam, _bai ->
                
                def normal_ids = normal_patient_ids.val
                paired: normal_ids.contains(meta.patient_id)
                unpaired: true
            }
            .unpaired
            .map { 
                _meta, _fastq_1, _fastq_2, bam, bai -> 
                def meta = [
                    id: _meta.sample_id,
                    patient_id: _meta.patient_id,
                    tumour_id: _meta.sample_id,
                    normal_id: null,
                    is_paired: false
                ]
                
                return [meta, bam, bai, null, null]
            }

        bam_tumour_only
            .count()
            .subscribe { count ->
                
                log.info("Total tumour-only samples = ${count}")
            }
        
        // Combine paired and unpaired samples
        bam_all_samples = bam_tumour_normal.mix(bam_tumour_only)
        
    emit:
        input_samples       = input_samples
        tumour_samples      = tumour_samples
        normal_samples      = normal_samples
        fastq               = fastq
        bam_tumour_normal   = bam_tumour_normal
        bam_tumour_only     = bam_tumour_only
        bam_all_samples     = bam_all_samples
}