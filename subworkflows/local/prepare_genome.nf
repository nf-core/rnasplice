//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_SUPPA_TPM        } from '../../modules/nf-core/modules/gunzip/main'

include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/modules/untar/main'

include { GFFREAD                           } from '../../modules/nf-core/modules/gffread/main'

include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/modules/star/genomegenerate/main'
include { STAR_GENOMEGENERATE_IGENOMES      } from '../../modules/local/star_genomegenerate_igenomes'

include { SALMON_INDEX                      } from '../../modules/nf-core/modules/salmon/index/main'

include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../modules/nf-core/modules/rsem/preparereference/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/modules/custom/getchromsizes/main'

include { GTF_GENE_FILTER                   } from '../../modules/local/gtf_gene_filter'

workflow PREPARE_GENOME {

    take:
    
    prepare_tool_indices // list   : tools to prepare indices for
    biotype              // string : if additional fasta file is provided biotype value to use when appending entries to GTF file
    is_aws_igenome       // boolean: whether the genome files are from AWS iGenomes

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)


    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], params.gff ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }


    //
    // Uncompress transcript fasta file / create if required
    //
    if (params.transcript_fasta) {
        if (params.transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], params.transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = file(params.transcript_fasta)
        } 
    } else {
        ch_filter_gtf       = GTF_GENE_FILTER ( ch_fasta, ch_gtf ).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_filter_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star' in prepare_tool_indices || 'star_salmon' in prepare_tool_indices) {
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], params.star_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = file(params.star_index)
            }
        } else {
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    if ('salmon' in prepare_tool_indices || 'star_salmon' in prepare_tool_indices) {
        if (params.salmon_index) {
            if (params.salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], params.salmon_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
            } else {
                ch_salmon_index = file(params.salmon_index)
            }
        } else {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    //
    // Gather Suppa tpm file 
    //
    ch_suppa_tpm = Channel.empty()
    if (params.suppa_tpm) {
        if (params.suppa_tpm.endsWith('.gz')) {
            ch_suppa_tpm = GUNZIP_SUPPA_TPM ( [ [:], params.suppa_tpm ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_SUPPA_TPM.out.versions)
        } else {
            ch_suppa_tpm = file(params.suppa_tpm)
        }
    }

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    fai              = ch_fai              //    path: genome.fai
    chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    gtf              = ch_gtf              //    path: genome.gtf
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    star_index       = ch_star_index       //    path: star/index/
    salmon_index     = ch_salmon_index     //    path: salmon/index/
    suppa_tpm        = ch_suppa_tpm        //    path: suppa.tpm

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
