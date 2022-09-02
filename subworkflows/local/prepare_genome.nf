//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/modules/gunzip/main'

include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/modules/untar/main'

include { GFFREAD                           } from '../../modules/nf-core/modules/gffread/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/modules/star/genomegenerate/main'
include { SALMON_INDEX                      } from '../../modules/nf-core/modules/salmon/index/main'

include { GTF2BED                      } from '../../modules/local/gtf2bed'
include { GTF_GENE_FILTER              } from '../../modules/local/gtf_gene_filter'

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
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], params.additional_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, biotype )
        ch_fasta    = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf      = CAT_ADDITIONAL_FASTA.out.gtf
        ch_versions = ch_versions.mix(CAT_ADDITIONAL_FASTA.out.versions)
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star_salmon' in prepare_tool_indices) {
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
    if ('salmon' in prepare_tool_indices) {
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

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    gtf              = ch_gtf              //    path: genome.gtf
    gene_bed         = ch_gene_bed         //    path: gene.bed
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    star_index       = ch_star_index       //    path: star/index/
    salmon_index     = ch_salmon_index     //    path: salmon/index/

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
