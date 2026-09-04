//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                          } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF                            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA               } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF_DEXSEQ                     } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_SUPPA_TPM                      } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_STAR_INDEX                       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX                     } from '../../../modules/nf-core/untar'

include { SAMTOOLS_FAIDX                                  } from '../../../modules/nf-core/samtools/faidx'
include { GFFREAD                                         } from '../../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE                             } from '../../../modules/nf-core/star/genomegenerate'
include { STAR_GENOMEGENERATE_IGENOMES                    } from '../../../modules/local/star_genomegenerate_igenomes'
include { SALMON_INDEX                                    } from '../../../modules/nf-core/salmon/index'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/rsem/preparereference'

include { GTFGENEFILTER                                   } from '../../../modules/local/gtfgenefilter'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE            } from '../../../modules/local/preprocess_transcripts_fasta_gencode'

workflow PREPARE_GENOME {
    take:
    fasta //      file: /path/to/genome.fasta
    gtf //      file: /path/to/genome.gtf
    gff //      file: /path/to/genome.gff
    transcript_fasta //      file: /path/to/transcript.fasta
    star_index // directory: /path/to/star/index/
    salmon_index // directory: /path/to/salmon/index/
    gff_dexseq //      file: /path/to/dexseq/genome.gff
    suppa_tpm //      file: /path/to/suppa/quant.tpm
    gencode //   boolean: whether gene annotation is from gencode
    is_aws_igenome //   boolean: whether the genome files are from AWS iGenomes

    main:

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA([[:], fasta])
        ch_fasta = GUNZIP_FASTA.out.gunzip
    }
    else {
        ch_fasta = channel.value([[:], file(fasta, checkIfExists: true)])
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            GUNZIP_GTF([[:], gtf])
            ch_gtf = GUNZIP_GTF.out.gunzip
        }
        else {
            ch_gtf = channel.value([[:], file(gtf, checkIfExists: true)])
        }
    }
    else if (gff) {
        if (gff.endsWith('.gz')) {
            GUNZIP_GFF([[:], gff])
            ch_gff = GUNZIP_GFF.out.gunzip
        }
        else {
            ch_gff = channel.value([[:], file(gff, checkIfExists: true)])
        }
        ch_gtf = GFFREAD(ch_gff, null).gtf
    }

    //
    // Uncompress transcript fasta file / create if required
    //
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            GUNZIP_TRANSCRIPT_FASTA([[:], transcript_fasta])
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA.out.gunzip
        }
        else {
            ch_transcript_fasta = channel.value([[:], file(transcript_fasta, checkIfExists: true)])
        }
        if (gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta)
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
        }
    }
    else {
        ch_filter_gtf = GTFGENEFILTER(ch_fasta, ch_gtf).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA(
            ch_fasta.map { _meta, fa -> fa },
            ch_filter_gtf.map { _meta, filter_gtf -> filter_gtf },
        ).transcript_fasta.map { fa -> [[:], fa] }
    }

    //
    // Create chromosome sizes file
    //
    SAMTOOLS_FAIDX(ch_fasta.map { meta, fa -> [meta, fa, []] }, true)
    ch_fai = SAMTOOLS_FAIDX.out.fai
    ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = channel.empty()
    if (params.source == 'fastq' && !params.skip_alignment && (params.aligner == 'star' || params.aligner == 'star_salmon')) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX([[:], star_index]).untar
            }
            else {
                ch_star_index = channel.value([[:], file(star_index, checkIfExists: true)])
            }
        }
        else {
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES(ch_fasta, ch_gtf).index
            }
            else {
                ch_star_index = STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
            }
        }
    }

    //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = channel.empty()
    if (params.source == 'fastq' && (params.pseudo_aligner == 'salmon' || params.aligner == 'star_salmon')) {
        if (salmon_index) {
            if (salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX([[:], salmon_index]).untar.map { _meta, index -> index }
            }
            else {
                ch_salmon_index = channel.value([[:], file(salmon_index, checkIfExists: true)])
            }
        }
        else {
            if (params.pseudo_aligner == 'salmon') {
                SALMON_INDEX(
                    ch_fasta.map { _meta, fa -> fa },
                    ch_transcript_fasta.map { _meta, tr -> tr },
                )
                ch_salmon_index = SALMON_INDEX.out.index
            }
        }
    }

    //
    // Uncompress DEXSeq GFF annotation file if required
    //
    ch_dexseq_gff = channel.empty()
    if (gff_dexseq) {
        if (gff_dexseq.endsWith('.gz')) {
            GUNZIP_GFF_DEXSEQ([[:], gff_dexseq])
            ch_dexseq_gff = GUNZIP_GFF_DEXSEQ.out.gunzip
        }
        else {
            ch_dexseq_gff = channel.value([[:], file(gff_dexseq, checkIfExists: true)])
        }
    }

    //
    // Uncompress SUPPA TPM file if required
    //
    ch_suppa_tpm = channel.empty()
    if (suppa_tpm) {
        if (suppa_tpm.endsWith('.gz')) {
            GUNZIP_SUPPA_TPM([[:], suppa_tpm])
            ch_suppa_tpm = GUNZIP_SUPPA_TPM.out.gunzip.map { _meta, tpm -> tpm }
        }
        else {
            ch_suppa_tpm = channel.value([[:], file(suppa_tpm, checkIfExists: true)])
        }
    }

    emit:
    fasta            = ch_fasta.map { _meta, fa -> fa } //    path: genome.fasta
    fai              = ch_fai //    path: genome.fai
    chrom_sizes      = ch_chrom_sizes //    path: genome.sizes
    gtf              = ch_gtf.map { _meta, out_gtf -> out_gtf } //    path: genome.gtf
    transcript_fasta = ch_transcript_fasta.map { _meta, fa -> fa } //    path: transcript.fasta
    star_index       = ch_star_index.map { _meta, index -> index } //    path: star/index/
    salmon_index     = ch_salmon_index //    path: salmon/index/
    dexseq_gff       = ch_dexseq_gff.map { _meta, dexseq_gff -> dexseq_gff } //    path: dexseq.gff
    suppa_tpm        = ch_suppa_tpm //    path: suppa.tpm
}
