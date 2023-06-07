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

include { GTF_GENE_FILTER                      } from '../../modules/local/gtf_gene_filter'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../modules/local/preprocess_transcripts_fasta_gencode.nf'

workflow PREPARE_GENOME {

    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf
    gff                  //      file: /path/to/genome.gff
    transcript_fasta     //      file: /path/to/transcript.fasta
    star_index           // directory: /path/to/star/index/
    salmon_index         // directory: /path/to/salmon/index/
    gff_dexseq           //      file: /path/to/dexseq/genome.gff
    suppa_tpm            //      file: /path/to/suppa/quant.tpm
    is_aws_igenome       //   boolean: whether the genome files are from AWS iGenomes
    prepare_tool_indices //      list: tools to prepare indices for
    step                 //     value: workflow step
    gencode              //   boolean: whether gene annotation is from gencode

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value(file(gff))
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress transcript fasta file / create if required
    //

    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_filter_gtf = GTF_GENE_FILTER ( ch_fasta, ch_gtf ).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_filter_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (('star' in prepare_tool_indices || 'star_salmon' in prepare_tool_indices) && (step == 'fastq')) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(star_index))
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
    if (('salmon' in prepare_tool_indices || 'star_salmon' in prepare_tool_indices) && (step == 'fastq')) {
        if (salmon_index) {
            if (salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
            } else {
                ch_salmon_index = Channel.value(file(salmon_index))
            }
        } else {
            if ('salmon' in prepare_tool_indices) {
                ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
                ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
            }
        }
    }

    //
    // Uncompress DEXSeq GFF annotation file if required
    //
    ch_dexseq_gff = Channel.empty()
    if (gff_dexseq) {
        if (gff_dexseq.endsWith('.gz')) {
            ch_dexseq_gff = GUNZIP_GFF_DEXSEQ ( [ [:], gff_dexseq ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF_DEXSEQ.out.versions)
        } else {
            ch_dexseq_gff = Channel.value(file(gff_dexseq))
        }
    }

    //
    // Uncompress SUPPA TPM file if required
    //
    ch_suppa_tpm = Channel.empty()
    if (suppa_tpm) {
        if (suppa_tpm.endsWith('.gz')) {
            ch_suppa_tpm = GUNZIP_SUPPA_TPM ( [ [:], suppa_tpm ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_SUPPA_TPM.out.versions)
        } else {
            ch_suppa_tpm = Channel.value(file(suppa_tpm))
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
