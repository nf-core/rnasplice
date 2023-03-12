//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF_DEXSEQ       } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_SUPPA_TPM        } from '../../modules/nf-core/gunzip/main'

include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/untar/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes/main'
include { GFFREAD                           } from '../../modules/nf-core/gffread/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_GENOMEGENERATE_IGENOMES      } from '../../modules/local/star_genomegenerate_igenomes'
include { SALMON_INDEX                      } from '../../modules/nf-core/salmon/index/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../modules/nf-core/rsem/preparereference/main'

include { GTF_GENE_FILTER                   } from '../../modules/local/gtf_gene_filter'

workflow PREPARE_GENOME {

    take:

    prepare_tool_indices // list   : tools to prepare indices for
    is_aws_igenome       // boolean: whether the genome files are from AWS iGenomes

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA ( [ [:], params.fasta ] )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.fromPath(params.fasta)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( [ [:], ch_fasta ] )
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            GUNZIP_GTF ( [ [:], params.gtf ] )
            ch_gtf      = GUNZIP_GTF.outgunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.fromPath(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF ( [ [:], params.gff ] )
            ch_gff      = GUNZIP_GFF.out.gunzip.map{ it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.fromPath(params.gff)
        }
        GFFREAD ( ch_gff )
        ch_gtf      = GFFREAD.out.gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (params.transcript_fasta) {
        if (params.transcript_fasta.endsWith('.gz')) {
            GUNZIP_TRANSCRIPT_FASTA ( [ [:], params.transcript_fasta ] )
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA.out.gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.fromPath(params.transcript_fasta)
        }
    } else {
        GTF_GENE_FILTER ( ch_fasta, ch_gtf )
        ch_filter_gtf       = GTF_GENE_FILTER.out.gtf
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_filter_gtf )
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA.out.transcript_fasta
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    // //
    // // Create chromosome sizes file
    // //
    // CUSTOM_GETCHROMSIZES ( [ [:], ch_fasta ] )
    // ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    // ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    // ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

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
                ch_star_index = Channel.fromPath(params.star_index)
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
                UNTAR_SALMON_INDEX ( [ [:], params.salmon_index ] )
                ch_salmon_index = UNTAR_SALMON_INDEX.out.untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
            } else {
                ch_salmon_index = Channel.fromPath(params.salmon_index)
            }
        } else {
            SALMON_INDEX ( ch_fasta, ch_transcript_fasta )
            ch_salmon_index = SALMON_INDEX.out.index
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    //
    // Uncompress DEXSeq GFF annotation file if required
    //
    ch_dexseq_gff = Channel.empty()
    if (params.gff_dexseq) {
        if (params.gff_dexseq.endsWith('.gz')) {
            GUNZIP_GFF_DEXSEQ ( [ [:], params.gff_dexseq ] )
            ch_dexseq_gff = GUNZIP_GFF_DEXSEQ.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF_DEXSEQ.out.versions)
        } else {
            ch_dexseq_gff = Channel.fromPath(params.gff_dexseq)
        }
    }

    //
    // Uncompress SUPPA TPM file if required
    //
    ch_suppa_tpm = Channel.empty()
    if (params.suppa_tpm) {
        if (params.suppa_tpm.endsWith('.gz')) {
            GUNZIP_SUPPA_TPM ( [ [:], params.suppa_tpm ] )
            ch_suppa_tpm = GUNZIP_SUPPA_TPM.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_SUPPA_TPM.out.versions)
        } else {
            ch_suppa_tpm = Channel.fromPath(params.suppa_tpm)
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
    dexseq_gff       = ch_dexseq_gff       //    path: dexseq.gff
    suppa_tpm        = ch_suppa_tpm        //    path: suppa.tpm

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
