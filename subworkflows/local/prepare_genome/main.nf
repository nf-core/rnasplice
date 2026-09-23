//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF_DEXSEQ } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_SUPPA_TPM } from '../../../modules/nf-core/gunzip'

include { UNTAR as UNTAR_STAR_INDEX } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_SALMON_INDEX } from '../../../modules/nf-core/untar'

include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx'
include { GFFREAD } from '../../../modules/nf-core/gffread'
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate'
include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/rsem/preparereference'

include { GTFGENEFILTER } from '../../../modules/local/gtfgenefilter'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../../modules/local/preprocess_transcripts_fasta_gencode'
include { STAR_GENOMEPARAMS_UPGRADE } from '../../../modules/local/star_genomeparams_upgrade'

workflow PREPARE_GENOME {
    take:
    fasta // string: path to the genome FASTA, may be gzipped
    gtf // string: path to the GTF annotation, may be gzipped
    gff // string: path to the GFF3 annotation, used when no GTF is given, may be gzipped
    transcript_fasta // string: path to the transcript FASTA, may be gzipped
    star_index // string: path to the STAR index directory, may be a .tar.gz archive
    salmon_index // string: path to the Salmon index directory, may be a .tar.gz archive
    gff_dexseq // string: path to the flattened DEXSeq GFF annotation, may be gzipped
    suppa_tpm // string: path to the SUPPA transcript TPM table, may be gzipped
    gencode // boolean: whether the annotation is from GENCODE
    source // string: type of input data [fastq, genome_bam, transcriptome_bam, salmon_results]
    aligner // string: genome aligner [star, star_salmon]
    pseudo_aligner // string: pseudo aligner [salmon]
    skip_alignment // boolean: whether the genome alignment is skipped

    main:

    //
    // MODULE: GUNZIP_FASTA
    //

    if (fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA([[:], file(fasta, checkIfExists: true)]).gunzip
    }
    else {
        ch_fasta = channel.value([[:], file(fasta, checkIfExists: true)])
    }

    //
    // MODULE: GUNZIP_GTF, GUNZIP_GFF and GFFREAD
    //

    // A GTF annotation is used as is, a GFF3 one is converted to GTF
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF([[:], file(gtf, checkIfExists: true)]).gunzip
        }
        else {
            ch_gtf = channel.value([[:], file(gtf, checkIfExists: true)])
        }
    }
    else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF([[:], file(gff, checkIfExists: true)]).gunzip
        }
        else {
            ch_gff = channel.value([[:], file(gff, checkIfExists: true)])
        }
        // GFFREAD names its output after `meta.id`, so give it the annotation name
        ch_gtf = GFFREAD(ch_gff.map { _meta, gff_file -> [[id: gff_file.baseName], gff_file] }, []).gtf
    }

    //
    // MODULE: GUNZIP_TRANSCRIPT_FASTA, PREPROCESS_TRANSCRIPTS_FASTA_GENCODE, GTFGENEFILTER and MAKE_TRANSCRIPTS_FASTA
    //

    // Without a transcript FASTA, one is built from the genome and the annotation, keeping
    // only the genes on sequences present in the genome FASTA
    if (transcript_fasta) {
        if (transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA([[:], file(transcript_fasta, checkIfExists: true)]).gunzip
        }
        else {
            ch_transcript_fasta = channel.value([[:], file(transcript_fasta, checkIfExists: true)])
        }
        if (gencode) {
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE(ch_transcript_fasta).fasta
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
    // MODULE: SAMTOOLS_FAIDX
    //

    SAMTOOLS_FAIDX(ch_fasta.map { meta, fa -> [meta, fa, []] }, true)

    //
    // MODULE: UNTAR_STAR_INDEX, STAR_GENOMEPARAMS_UPGRADE and STAR_GENOMEGENERATE
    //

    ch_star_index = channel.empty()
    if (source == 'fastq' && !skip_alignment && aligner in ['star', 'star_salmon']) {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index_raw = UNTAR_STAR_INDEX([[:], file(star_index, checkIfExists: true)]).untar
            }
            else {
                ch_star_index_raw = channel.value([[:], file(star_index, checkIfExists: true)])
            }

            // A supplied index may have been built with STAR 2.6.x, as the AWS iGenomes ones were,
            // which STAR 2.7.4a and later refuse to read. `STAR_GENOMEPARAMS_UPGRADE` rewrites the
            // `genomeParameters.txt` metadata that changed and leaves a modern index untouched.
            ch_star_index = STAR_GENOMEPARAMS_UPGRADE(ch_star_index_raw).index
        }
        else {
            ch_star_index = STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
        }
    }

    //
    // MODULE: UNTAR_SALMON_INDEX and SALMON_INDEX
    //

    // `star_salmon` quantifies the STAR transcriptome alignments, which needs no index, so an
    // index is only built for the `salmon` pseudo aligner
    ch_salmon_index = channel.empty()
    if (source == 'fastq' && (pseudo_aligner == 'salmon' || aligner == 'star_salmon')) {
        if (salmon_index) {
            if (salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR_SALMON_INDEX([[:], file(salmon_index, checkIfExists: true)]).untar
            }
            else {
                ch_salmon_index = channel.value([[:], file(salmon_index, checkIfExists: true)])
            }
        }
        else if (pseudo_aligner == 'salmon') {
            ch_salmon_index = SALMON_INDEX(
                ch_fasta.map { _meta, fa -> fa },
                ch_transcript_fasta.map { _meta, tr -> tr },
            ).index.map { index -> [[:], index] }
        }
    }

    //
    // MODULE: GUNZIP_GFF_DEXSEQ
    //

    ch_dexseq_gff = channel.empty()
    if (gff_dexseq) {
        if (gff_dexseq.endsWith('.gz')) {
            ch_dexseq_gff = GUNZIP_GFF_DEXSEQ([[:], file(gff_dexseq, checkIfExists: true)]).gunzip
        }
        else {
            ch_dexseq_gff = channel.value([[:], file(gff_dexseq, checkIfExists: true)])
        }
    }

    //
    // MODULE: GUNZIP_SUPPA_TPM
    //

    ch_suppa_tpm = channel.empty()
    if (suppa_tpm) {
        if (suppa_tpm.endsWith('.gz')) {
            ch_suppa_tpm = GUNZIP_SUPPA_TPM([[:], file(suppa_tpm, checkIfExists: true)]).gunzip
        }
        else {
            ch_suppa_tpm = channel.value([[:], file(suppa_tpm, checkIfExists: true)])
        }
    }

    emit:
    fasta = ch_fasta.map { _meta, fa -> fa } // channel: path(genome.fasta)
    fai = SAMTOOLS_FAIDX.out.fai // channel: [ val(meta), path(genome.fasta.fai) ]
    chrom_sizes = SAMTOOLS_FAIDX.out.sizes // channel: [ val(meta), path(genome.fasta.sizes) ]
    gtf = ch_gtf.map { _meta, out_gtf -> out_gtf } // channel: path(genome.gtf)
    transcript_fasta = ch_transcript_fasta.map { _meta, fa -> fa } // channel: path(transcripts.fasta)
    star_index = ch_star_index.map { _meta, index -> index } // channel: path(star/index/)
    salmon_index = ch_salmon_index.map { _meta, index -> index } // channel: path(salmon/index/)
    dexseq_gff = ch_dexseq_gff.map { _meta, dexseq_gff -> dexseq_gff } // channel: path(dexseq.gff)
    suppa_tpm = ch_suppa_tpm.map { _meta, tpm -> tpm } // channel: path(suppa.tpm)
}
