//
// Differential splicing with SUPPA, per local event and per isoform
//

include { SUPPA_GENERATEEVENTS as SUPPA_GENERATEEVENTS_IOE } from '../../../modules/nf-core/suppa/generateevents'
include { SUPPA_GENERATEEVENTS as SUPPA_GENERATEEVENTS_IOI } from '../../../modules/nf-core/suppa/generateevents'
include { SUPPA_PSIPEREVENT } from '../../../modules/nf-core/suppa/psiperevent'
include { SUPPA_PSIPERISOFORM } from '../../../modules/nf-core/suppa/psiperisoform'
include { SUPPA_DIFFSPLICE } from '../../../modules/nf-core/suppa/diffsplice'
include { SUPPA_CLUSTEREVENTS } from '../../../modules/nf-core/suppa/clusterevents'
include { MERGEEVENTS } from '../../../modules/local/mergeevents'
include { SPLIT_FILES as SPLIT_FILES_TPM } from '../../../modules/local/splitfiles'
include { SPLIT_FILES as SPLIT_FILES_PSI } from '../../../modules/local/splitfiles'
include { CLUSTERGROUPS } from '../../../modules/local/clustergroups'

workflow SUPPA {
    take:
    ch_gtf // channel: [ val(meta), path(gtf) ]
    ch_tpm // channel: [ val(meta), path(tpm) ]
    ch_samples // channel: [ val(sample_id), val(condition) ], in samplesheet order
    ch_contrastsheet // channel: [ contrast: val(contrast), treatment: val(condition), control: val(condition) ]
    suppa_per_local_event // boolean: params.suppa_per_local_event
    generateevents_boundary // string: params.generateevents_boundary
    generateevents_threshold // integer: params.generateevents_threshold
    generateevents_exon_length // integer: params.generateevents_exon_length
    generateevents_event_type // string: params.generateevents_event_type
    generateevents_pool_genes // boolean: params.generateevents_pool_genes
    psiperevent_total_filter // integer: params.psiperevent_total_filter
    diffsplice_local_event // boolean: params.diffsplice_local_event
    diffsplice_transcript_event // boolean: params.diffsplice_isoform
    diffsplice_method // string: params.diffsplice_method
    diffsplice_area // integer: params.diffsplice_area
    diffsplice_lower_bound // float: params.diffsplice_lower_bound
    diffsplice_alpha // float: params.diffsplice_alpha
    diffsplice_tpm_threshold // float: params.diffsplice_tpm_threshold
    diffsplice_nan_threshold // float: params.diffsplice_nan_threshold
    diffsplice_gene_correction // boolean: params.diffsplice_gene_correction
    diffsplice_paired // boolean: params.diffsplice_paired
    diffsplice_median // boolean: params.diffsplice_median
    clusterevents_local_event // boolean: params.clusterevents_local_event
    clusterevents_transcript_event // boolean: params.clusterevents_isoform
    clusterevents_dpsithreshold // float: params.clusterevents_dpsithreshold
    clusterevents_eps // float: params.clusterevents_eps
    clusterevents_metric // string: params.clusterevents_metric
    clusterevents_min_pts // integer: params.clusterevents_min_pts
    clusterevents_method // string: params.clusterevents_method
    clusterevents_sigthreshold // float: params.clusterevents_sigthreshold
    clusterevents_separation // float: params.clusterevents_separation
    suppa_per_isoform // boolean: params.suppa_per_isoform

    main:

    // Sample ids of each condition, in samplesheet order: [ condition, [ sample_ids ] ].
    // The TPM and PSI files of a condition get their columns in this order, which is
    // what pairs the samples of two conditions in the paired diffSplice model
    ch_condition_samples = ch_samples
        .map { sample_id, condition -> [condition, sample_id] }
        .unique()
        .groupTuple()

    //
    // MODULE: Split the TPM file into one file per condition
    //
    SPLIT_FILES_TPM(
        ch_tpm.combine(ch_condition_samples).map { meta, tpm, condition, sample_ids ->
            [meta + [id: "${meta.id}_${condition}".toString(), condition: condition], tpm, sample_ids]
        },
        'tpm',
    )

    ch_events_ioe = channel.empty()
    ch_events_psi = channel.empty()

    if (suppa_per_local_event) {

        //
        // MODULE: Local events of the annotation, one file per event type
        //
        SUPPA_GENERATEEVENTS_IOE(
            ch_gtf,
            'ioe',
            generateevents_pool_genes,
            generateevents_event_type,
            generateevents_boundary,
            generateevents_threshold,
            generateevents_exon_length,
        )

        //
        // MODULE: Merge the event types into one file
        //
        MERGEEVENTS(SUPPA_GENERATEEVENTS_IOE.out.events)
        ch_events_ioe = MERGEEVENTS.out.ioe

        //
        // MODULE: PSI value of each local event
        //
        SUPPA_PSIPEREVENT(
            ch_tpm,
            ch_events_ioe,
            psiperevent_total_filter,
        )
        ch_events_psi = SUPPA_PSIPEREVENT.out.psi
    }

    ch_events_ioi = channel.empty()
    ch_isoform_psi = channel.empty()

    if (suppa_per_isoform) {

        //
        // MODULE: Transcript events of the annotation
        //
        SUPPA_GENERATEEVENTS_IOI(
            ch_gtf,
            'ioi',
            generateevents_pool_genes,
            generateevents_event_type,
            generateevents_boundary,
            generateevents_threshold,
            generateevents_exon_length,
        )
        ch_events_ioi = SUPPA_GENERATEEVENTS_IOI.out.events

        //
        // MODULE: PSI value of each isoform
        //
        SUPPA_PSIPERISOFORM(ch_tpm, ch_gtf)
        ch_isoform_psi = SUPPA_PSIPERISOFORM.out.psi
    }

    // From here on, both levels go through the same steps. The level, local or
    // transcript, is in the meta map and starts the name of every file:
    // [ level, events ] and [ level, psi ]
    ch_level_events = ch_events_ioe
        .map { _meta, ioe -> ['local', ioe] }
        .mix(ch_events_ioi.map { _meta, ioi -> ['transcript', ioi] })

    ch_level_psi = ch_events_psi
        .map { _meta, psi -> ['local', psi] }
        .mix(ch_isoform_psi.map { _meta, psi -> ['transcript', psi] })

    //
    // MODULE: Split the PSI files into one file per condition
    //
    SPLIT_FILES_PSI(
        ch_level_psi.combine(ch_condition_samples).map { level, psi, condition, sample_ids ->
            [[id: "${level}_${condition}".toString(), level: level, condition: condition], psi, sample_ids]
        },
        'psi',
    )

    // TPM and PSI files of each condition: [ [ level, condition ], tpm, psi ]
    ch_condition_tpm_psi = SPLIT_FILES_PSI.out.split
        .map { meta, psi -> [meta.condition, meta.level, psi] }
        .combine(SPLIT_FILES_TPM.out.split.map { meta, tpm -> [meta.condition, tpm] }, by: 0)
        .map { condition, level, psi, tpm -> [[level, condition], tpm, psi] }

    // The events, and the files of the treatment and of the control, of each contrast:
    // [ meta, events, tpm1, psi1, tpm2, psi2 ]
    ch_diffsplice = ch_level_events
        .filter { level, _events -> level == 'local' ? diffsplice_local_event : diffsplice_transcript_event }
        .combine(ch_contrastsheet)
        .map { level, events, contrast ->
            def meta = [
                id: "${level}_${contrast.treatment}-${contrast.control}".toString(),
                level: level,
                treatment: contrast.treatment,
                control: contrast.control,
            ]
            [[level, contrast.treatment], meta, events]
        }
        .combine(ch_condition_tpm_psi, by: 0)
        .map { _key, meta, events, tpm1, psi1 -> [[meta.level, meta.control], meta, events, tpm1, psi1] }
        .combine(ch_condition_tpm_psi, by: 0)
        .map { _key, meta, events, tpm1, psi1, tpm2, psi2 -> [meta, events, tpm1, psi1, tpm2, psi2] }

    //
    // MODULE: Differential splicing between the conditions of each contrast
    //
    SUPPA_DIFFSPLICE(
        ch_diffsplice.map { meta, events, _tpm1, _psi1, _tpm2, _psi2 -> [meta, events] },
        ch_diffsplice.map { meta, _events, tpm1, psi1, _tpm2, _psi2 -> [meta, meta.treatment, tpm1, psi1] },
        ch_diffsplice.map { meta, _events, _tpm1, _psi1, tpm2, psi2 -> [meta, meta.control, tpm2, psi2] },
        diffsplice_method,
        diffsplice_area,
        diffsplice_lower_bound,
        diffsplice_paired,
        diffsplice_gene_correction,
        diffsplice_alpha,
        false,
        false,
        diffsplice_median,
        diffsplice_tpm_threshold,
        diffsplice_nan_threshold,
    )

    ch_cluster_psivec = SUPPA_DIFFSPLICE.out.psivec.filter { meta, _psivec -> meta.level == 'local' ? clusterevents_local_event : clusterevents_transcript_event }

    //
    // MODULE: Column ranges of the two conditions in the PSI vector file
    //
    CLUSTERGROUPS(ch_cluster_psivec)

    // [ meta, dpsi, psivec, ranges ]
    ch_clusterevents = SUPPA_DIFFSPLICE.out.dpsi
        .join(ch_cluster_psivec)
        .join(CLUSTERGROUPS.out.ranges)

    //
    // MODULE: Cluster the events by their PSI values across the samples
    //
    SUPPA_CLUSTEREVENTS(
        ch_clusterevents.map { meta, dpsi, psivec, _ranges -> [meta, dpsi, psivec] },
        clusterevents_sigthreshold,
        clusterevents_dpsithreshold,
        clusterevents_eps,
        clusterevents_metric,
        clusterevents_separation,
        clusterevents_min_pts,
        ch_clusterevents.map { _meta, _dpsi, _psivec, ranges -> ranges },
        clusterevents_method,
    )

    emit:
    ioe_events = ch_events_ioe // channel: [ val(meta), path(ioe) ]
    ioi_events = ch_events_ioi // channel: [ val(meta), path(ioi) ]
    suppa_local_psi = ch_events_psi // channel: [ val(meta), path(psi) ]
    suppa_isoform_psi = ch_isoform_psi // channel: [ val(meta), path(psi) ]
    split_suppa_tpms = SPLIT_FILES_TPM.out.split // channel: [ val(meta), path(tpm) ], one per condition
    split_suppa_local_psi = SPLIT_FILES_PSI.out.split.filter { meta, _psi -> meta.level == 'local' } // channel: [ val(meta), path(psi) ], one per condition
    split_suppa_isoform_psi = SPLIT_FILES_PSI.out.split.filter { meta, _psi -> meta.level == 'transcript' } // channel: [ val(meta), path(psi) ], one per condition
    dpsi_local = SUPPA_DIFFSPLICE.out.dpsi.filter { meta, _dpsi -> meta.level == 'local' } // channel: [ val(meta), path(dpsi) ]
    psivec_local = SUPPA_DIFFSPLICE.out.psivec.filter { meta, _psivec -> meta.level == 'local' } // channel: [ val(meta), path(psivec) ]
    dpsi_isoform = SUPPA_DIFFSPLICE.out.dpsi.filter { meta, _dpsi -> meta.level == 'transcript' } // channel: [ val(meta), path(dpsi) ]
    psivec_isoform = SUPPA_DIFFSPLICE.out.psivec.filter { meta, _psivec -> meta.level == 'transcript' } // channel: [ val(meta), path(psivec) ]
    groups = CLUSTERGROUPS.out.groups // channel: [ val(meta), path(txt) ]
    cluster_vec_local = SUPPA_CLUSTEREVENTS.out.clustvec.filter { meta, _clustvec -> meta.level == 'local' } // channel: [ val(meta), path(clustvec) ]
    cluster_vec_isoform = SUPPA_CLUSTEREVENTS.out.clustvec.filter { meta, _clustvec -> meta.level == 'transcript' } // channel: [ val(meta), path(clustvec) ]
}
