
include { VARIANT_ANNOTATION } from '../../subworkflows/mutation_calling/variant_annotation.nf'

params.vcf     = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_filter/DFSP-001-T/DFSP-001-T.final.vcf.gz"
params.vcf_tbi = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_filter/DFSP-001-T/DFSP-001-T.final.vcf.gz.tbi"

params.fasta                    = params.genomes[params.genome]?.fasta
params.fai                      = params.genomes[params.genome]?.fai
params.dict                     = params.genomes[params.genome]?.dict
params.dbsnp                    = params.genomes[params.genome]?.dbsnp
params.dbsnp_tbi                = params.genomes[params.genome]?.dbsnp_tbi
params.germline_resource        = params.genomes[params.genome]?.germline_resource
params.germline_resource_tbi    = params.genomes[params.genome]?.germline_resource_tbi
params.panel_of_normals         = params.genomes[params.genome]?.pon
params.panel_of_normals_tbi     = params.genomes[params.genome]?.pon_tbi
params.pileup_variants          = params.genomes[params.genome]?.contamination_variants
params.pileup_variants_tbi      = params.genomes[params.genome]?.contamination_variants_tbi
params.intervals                = params.genomes[params.genome]?.intervals
params.bait_intervals           = params.genomes[params.genome]?.bait_intervals
params.target_intervals         = params.genomes[params.genome]?.target_intervals
params.targets                  = params.genomes[params.genome]?.targets

params.funcotator_resources     = params.genomes[params.genome]?.funcotator_resources
params.funcotator_ref_version   =params.genomes[params.genome]?.funcotator_ref_version

params.annovar_db               = params.genomes[params.genome]?.annovar_db
params.annovar_buildver         = params.genomes[params.genome]?.annovar_buildver
params.annovar_protocol         = params.genomes[params.genome]?.annovar_protocol
params.annovar_operation        = params.genomes[params.genome]?.annovar_operation
params.annovar_xreffile         = params.genomes[params.genome]?.annovar_xreffile
params.repeatmasker             = params.genomes[params.genome]?.repeatmasker
params.blacklist                = params.genomes[params.genome]?.blacklist

workflow {

    fasta                   = params.fasta ? file(params.fasta) : null
    fai                     = params.fai ? file(params.fai) : null
    dict                    = params.dict ? file(params.dict) : null
    dbsnp                   = params.dbsnp ? file(params.dbsnp) : null
    dbsnp_tbi               = params.dbsnp_tbi ? file(params.dbsnp_tbi) : null
    germline_resource       = params.germline_resource ? file(params.germline_resource) : null
    germline_resource_tbi   = params.germline_resource_tbi ? file(params.germline_resource_tbi) : null
    panel_of_normals        = params.panel_of_normals ? file(params.panel_of_normals) : null
    panel_of_normals_tbi    = params.panel_of_normals_tbi ? file(params.panel_of_normals_tbi) : null
    pileup_variants         = params.pileup_variants ? file(params.pileup_variants) : null
    pileup_variants_tbi     = params.pileup_variants_tbi ? file(params.pileup_variants_tbi) : null
    intervals               = params.intervals ? file(params.intervals) : null
    bait_intervals          = params.bait_intervals ? file(params.bait_intervals) : null
    target_intervals        = params.target_intervals ? file(params.target_intervals) : null
    targets                 = params.targets ? file(params.targets) : null
    repeatmasker            = params.repeatmasker ? file(params.repeatmasker) : null
    blacklist               = params.blacklist ? file(params.blacklist) : null

    funcotator_resources    = params.funcotator_resources ? file(params.funcotator_resources) : null
    annovar_db              = params.annovar_db ? file(params.annovar_db) : null
    annovar_buildver        = params.annovar_buildver
    annovar_protocol        = params.annovar_protocol
    annovar_operation       = params.annovar_operation
    annovar_xreffile        = params.annovar_xreffile ? file(params.annovar_xreffile) : null

    vcf = file(params.vcf)
    vcf_tbi = file(params.vcf_tbi)
    funcotator_ref_version = params.funcotator_ref_version

    def meta = [
        id: 'DFSP-001-T',
        patient_id: 'DFSP-001',
        sample_id: 'DFSP-001-T',
        normal_id: 'DFSP-001-N',
        is_paired: true
    ]

    input_vcf = Channel.of([meta, vcf, vcf_tbi])
    input_vcf.view()

    VARIANT_ANNOTATION(
        input_vcf,
        fasta,
        fai,
        dict,
        intervals,
        funcotator_resources,
        funcotator_ref_version,
        annovar_db,
        annovar_buildver,
        annovar_protocol,
        annovar_operation,
        annovar_xreffile
    )
}