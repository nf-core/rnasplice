//
// Tx2gene and Tximport subworkflow
//

include { GFFREAD_TX2GENE  } from '../../modules/local/gffread_tx2gene'
include { TXIMPORT         } from '../../modules/local/tximport'
include { UNTAR            } from '../../modules/nf-core/untar'

workflow TX2GENE_TXIMPORT {

    take:

    salmon_results   // tuple [meta, salmon path] when starting from salmon or path: "salmon/*" (SALMON_QUANT.out.results.collect{it[1]}) other cases
    gtf              // channel: /path/to/genome.gtf

    main:

    ch_versions = Channel.empty()

    //
    // Extract archives (if necessary)
    //

    salmon_results
    .map {
        meta, prefix ->
            tgz = prefix[0].toString().endsWith(".tar.gz") ? true : false
            [ meta + [tgz: tgz], prefix ]
    }
    .branch{
        tar: it[0].tgz == true
        dir: it[0].tgz == false
    }
    .set{ salmon_results }
    UNTAR ( salmon_results.tar )
    salmon_results = salmon_results.dir.mix(UNTAR.out.untar)

    //
    // Quantify and merge counts across samples
    //

    GFFREAD_TX2GENE ( gtf )

    ch_versions = ch_versions.mix(GFFREAD_TX2GENE.out.versions)

    tx2gene = GFFREAD_TX2GENE.out.tx2gene

    TXIMPORT ( salmon_results.collect{it[1]}, tx2gene )

    ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    emit:

    tx2gene                                                                         // path: *.tx2gene.tsv

    txi                             = TXIMPORT.out.txi                              // path: txi.rds
    txi_s                           = TXIMPORT.out.txi_s                            // path: txi.s.rds
    txi_ls                          = TXIMPORT.out.txi_ls                           // path: txi.ls.rds
    txi_dtu                         = TXIMPORT.out.txi_dtu                          // path: txi.dtu.rds

    gi                              = TXIMPORT.out.gi                               // path: gi.rds
    gi_s                            = TXIMPORT.out.gi_s                             // path: gi.s.rds
    gi_ls                           = TXIMPORT.out.gi_ls                            // path: gi.ls.rds

    tpm_gene                        = TXIMPORT.out.tpm_gene                         // path: *.gene_tpm.tsv
    counts_gene                     = TXIMPORT.out.counts_gene                      // path: *.gene_counts.tsv
    tpm_gene_scaled                 = TXIMPORT.out.tpm_gene_scaled                  // path: *.gene_tpm_scaled.tsv
    counts_gene_scaled              = TXIMPORT.out.counts_gene_scaled               // path: *.gene_counts_scaled.tsv
    tpm_gene_length_scaled          = TXIMPORT.out.tpm_gene_length_scaled           // path: *.gene_tpm_length_scaled.tsv
    counts_gene_length_scaled       = TXIMPORT.out.counts_gene_length_scaled        // path: *.gene_counts_length_scaled.tsv

    tpm_transcript                  = TXIMPORT.out.tpm_transcript                   // path: *.transcript_tpm.tsv
    counts_transcript               = TXIMPORT.out.counts_transcript                // path: *.transcript_counts.tsv
    tpm_transcript_scaled           = TXIMPORT.out.tpm_transcript_scaled            // path: *.transcript_tpm_scaled.tsv
    counts_transcript_scaled        = TXIMPORT.out.counts_transcript_scaled         // path: *.transcript_counts_scaled.tsv
    tpm_transcript_length_scaled    = TXIMPORT.out.tpm_transcript_length_scaled     // path: *.transcript_tpm_length_scaled.tsv
    counts_transcript_length_scaled = TXIMPORT.out.counts_transcript_length_scaled  // path: *.transcript_counts_length_scaled.tsv
    tpm_transcript_dtu_scaled       = TXIMPORT.out.tpm_transcript_dtu_scaled        // path: *.transcript_tpm_dtu_scaled.tsv
    counts_transcript_dtu_scaled    = TXIMPORT.out.counts_transcript_dtu_scaled     // path: *.transcript_counts_dtu_scaled.tsv

    tximport_tx2gene                = TXIMPORT.out.tximport_tx2gene                 // path: tximport.tx2gene.tsv

    suppa_tpm                       = TXIMPORT.out.suppa_tpm                        // path: suppa_tpm.txt

    versions                        = ch_versions                                   // channel: [ versions.yml ]
}
