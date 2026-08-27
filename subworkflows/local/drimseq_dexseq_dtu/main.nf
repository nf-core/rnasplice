//
// Dexseq DTU subworkflow
//

include { DRIMSEQ_FILTER  } from '../../../modules/local/drimseq_filter'
include { DEXSEQ_DTU      } from '../../../modules/local/dexseq/dtu'
include { STAGER          } from '../../../modules/local/stager'

workflow DRIMSEQ_DEXSEQ_DTU {

    take:

    txi                     // path: *.txi*.rds (either txi.s.rds or txi.dtu.rds)
    tximport_tx2gene        // path: tximport.tx2gene.tsv
    samplesheet             // path: /path/to/samplesheet.csv
    contrastsheet           // path: contrastsheet
    min_samps_gene_expr     // params.min_samps_gene_expr
    min_samps_feature_expr  // params.min_samps_feature_expr
    min_samps_feature_prop  // params.min_samps_feature_prop
    min_feature_expr        // params.min_feature_expr
    min_feature_prop        // params.min_feature_prop
    min_gene_expr           // params.min_gene_expr

    main:

    //
    // DEXSEQ FILTER
    //

    DRIMSEQ_FILTER (
        txi,
        tximport_tx2gene,
        samplesheet,
        min_samps_gene_expr,
        min_samps_feature_expr,
        min_samps_feature_prop,
        min_feature_expr,
        min_feature_prop,
        min_gene_expr
    )

    //
    // DEXSEQ DTU
    //

    DEXSEQ_DTU (
        DRIMSEQ_FILTER.out.drimseq_samples_tsv,
        DRIMSEQ_FILTER.out.drimseq_counts_tsv,
        contrastsheet
    )

    //
    // Join feature and gene channels by contrast value (extracted from filename)
    //

    ch_dexseq_feature_tsv = DEXSEQ_DTU.out.dexseq_exon_results_tsv
        .flatten()
        .map { it ->
            [ [ id: it.baseName.toString().replaceAll("DEXSeqResults.", "") ], it ]
        }

    ch_dexseq_gene_tsv = DEXSEQ_DTU.out.dexseq_gene_results_tsv
        .flatten()
        .map { it ->
            [ [ id: it.baseName.toString().replaceAll("perGeneQValue.", "") ], it ]
        }

    ch_dexseq_feature_gene_tsv = ch_dexseq_feature_tsv.join(ch_dexseq_gene_tsv)

    //
    // STAGER
    //

    def analysis_type = 'dexseq'

    STAGER (
        ch_dexseq_feature_gene_tsv,
        analysis_type
    )

    emit:

    drimseq_dataset_rds      = DRIMSEQ_FILTER.out.drimseq_dataset_rds
    drimseq_samples_tsv      = DRIMSEQ_FILTER.out.drimseq_samples_tsv
    drimseq_counts_tsv       = DRIMSEQ_FILTER.out.drimseq_counts_tsv

    dexseq_exon_dataset_rds  = DEXSEQ_DTU.out.dexseq_exon_dataset_rds
    dexseq_exon_results_rds  = DEXSEQ_DTU.out.dexseq_exon_results_rds
    dexseq_exon_results_tsv  = DEXSEQ_DTU.out.dexseq_exon_results_tsv
    dexseq_gene_results_rds  = DEXSEQ_DTU.out.dexseq_gene_results_rds
    dexseq_gene_results_tsv  = DEXSEQ_DTU.out.dexseq_gene_results_tsv

    stager_rds               = STAGER.out.stager_rds
    stager_padj_rds          = STAGER.out.stager_padj_rds
    stager_padj_tsv          = STAGER.out.stager_padj_tsv
}
