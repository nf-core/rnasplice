//
// SUPPA Subworkflow
//

include { SUPPA_GENERATEEVENTS as SUPPA_GENERATEEVENTS_IOE } from '../../../modules/nf-core/suppa/generateevents'
include { SUPPA_GENERATEEVENTS as SUPPA_GENERATEEVENTS_IOI } from '../../../modules/nf-core/suppa/generateevents'

include { SUPPA_PSIPEREVENT                                } from '../../../modules/nf-core/suppa/psiperevent'
include { SUPPA_PSIPERISOFORM                              } from '../../../modules/nf-core/suppa/psiperisoform'

include { SPLIT_FILES as SPLIT_FILES_TPM                   } from '../../../modules/local/splitfiles'
include { SPLIT_FILES as SPLIT_FILES_IOE                   } from '../../../modules/local/splitfiles'
include { SPLIT_FILES as SPLIT_FILES_IOI                   } from '../../../modules/local/splitfiles'

include { SUPPA_DIFFSPLICE as DIFFSPLICE_IOE               } from '../../../modules/nf-core/suppa/diffsplice'
include { SUPPA_DIFFSPLICE as DIFFSPLICE_IOI               } from '../../../modules/nf-core/suppa/diffsplice'

include { CLUSTERGROUPS as CLUSTERGROUPS_IOE               } from '../../../modules/local/clustergroups'
include { CLUSTERGROUPS as CLUSTERGROUPS_IOI               } from '../../../modules/local/clustergroups'
include { MERGEEVENTS                                      } from '../../../modules/local/mergeevents'

include { SUPPA_CLUSTEREVENTS as CLUSTEREVENTS_IOE         } from '../../../modules/nf-core/suppa/clusterevents'
include { SUPPA_CLUSTEREVENTS as CLUSTEREVENTS_IOI         } from '../../../modules/nf-core/suppa/clusterevents'

workflow SUPPA {

    take:

    ch_gtf                         // [ meta, gtf ]
    ch_tpm                         // [ meta, tpm ]
    ch_samplesheet                 // path(samplesheet)
    ch_contrastsheet               // channel: [ contrast, treatment, control ]
    suppa_per_local_event          // params.suppa_per_local_event
    generateevents_boundary        // params.generateevents_boundary
    generateevents_threshold       // params.generateevents_threshold
    generateevents_exon_length     // params.generateevents_exon_length
    generateevents_event_type      // params.generateevents_event_type
    generateevents_pool_genes      // params.generateevents_pool_genes
    psiperevent_total_filter       // params.psiperevent_total_filter
    diffsplice_local_event         // params.diffsplice_local_event
    diffsplice_transcript_event    // params.diffsplice_isoform
    diffsplice_method              // params.diffsplice_method
    diffsplice_area                // params.diffsplice_area
    diffsplice_lower_bound         // params.diffsplice_lower_bound
    diffsplice_alpha               // params.diffsplice_alpha
    diffsplice_tpm_threshold       // params.diffsplice_tpm_threshold
    diffsplice_nan_threshold       // params.diffsplice_nan_threshold
    diffsplice_gene_correction     // params.diffsplice_gene_correction
    diffsplice_paired              // params.diffsplice_paired
    diffsplice_median              // params.diffsplice_median
    clusterevents_local_event      // params.clusterevents_local_event
    clusterevents_transcript_event // params.clusterevents_transcript_event
    clusterevents_dpsithreshold    // params.clusterevents_dpsithreshold
    clusterevents_eps              // params.clusterevents_eps
    clusterevents_metric           // params.clusterevents_metric
    clusterevents_min_pts          // params.clusterevents_min_pts
    clusterevents_method           // params.clusterevents_method
    clusterevents_sigthreshold     // params.clusterevents_sigthreshold
    clusterevents_separation       // params.clusterevents_separation
    suppa_per_isoform              // params.suppa_per_isoform

    main:

    // Split the tpm file (contains all samples) into individual files based on condition

    SPLIT_FILES_TPM (
        ch_tpm,
        ch_samplesheet,
        ".tpm",
        false
    )

    ch_split_tpms = SPLIT_FILES_TPM.out.tpms

    // If per AS local analysis:

    ch_events_ioe             = channel.empty()
    ch_events_psi             = channel.empty()
    ch_split_events_psi       = channel.empty()

    ch_events_dpsi            = channel.empty()
    ch_events_psivec          = channel.empty()

    ch_groups_ioe             = channel.empty()
    ch_events_clustvec        = channel.empty()

    if (suppa_per_local_event) {

        // Generate AS events on the GTF - Local events

        SUPPA_GENERATEEVENTS_IOE (
            ch_gtf,
            'ioe',
            generateevents_pool_genes,
            generateevents_event_type,
            generateevents_boundary,
            generateevents_threshold,
            generateevents_exon_length
        )

        ch_events_ioe  = MERGEEVENTS(SUPPA_GENERATEEVENTS_IOE.out.events).ioe

        // Calculate the psi values of Local events (using events file and TPM)

        SUPPA_PSIPEREVENT (
            ch_tpm,
            ch_events_ioe,
            psiperevent_total_filter
        )

        ch_events_psi = SUPPA_PSIPEREVENT.out.psi

        // Split the PSI files between the conditions

        SPLIT_FILES_IOE (
            ch_events_psi,
            ch_samplesheet,
            '.psi',
            true
        )

        ch_split_events_psi = SPLIT_FILES_IOE.out.psis

        // Calculate differential analysis between conditions

        if (diffsplice_local_event) {

            // Create contrasts channel

            ch_suppa_local_contrasts_raw = ch_contrastsheet

            // Add TPM files to contrasts channel

            SPLIT_FILES_TPM.out.tpms
                .flatten()
                .map { it -> [it.baseName, it ] }
                .multiMap { base, file ->
                    for_tpm1: [ base, file ]
                    for_tpm2: [ base, file ]
                }
                .set { ch_suppa_local_tpm_conditions }

            ch_suppa_local_contrasts_tpm1 = ch_suppa_local_contrasts_raw
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_local_tpm_conditions.for_tpm1, by: 0 )
                .map { it -> it[1] + ['tpm1': it[2]] }

            ch_suppa_local_contrasts_tpm2 = ch_suppa_local_contrasts_tpm1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_local_tpm_conditions.for_tpm2, by: 0 )
                .map { it -> it[1] + ['tpm2': it[2]] }

            // Add PSI files to contrasts channel

            SPLIT_FILES_IOE.out.psis
                .flatten()
                .map { it -> [ it.baseName.toString().replaceAll('local_', ''), it ] }
                .multiMap { base, file ->
                    for_psi1: [ base, file ]
                    for_psi2: [ base, file ]
                }
                .set { ch_suppa_local_psi_conditions }

            ch_suppa_local_contrasts_psi1 = ch_suppa_local_contrasts_tpm2
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_local_psi_conditions.for_psi1, by: 0 )
                .map { it -> it[1] + ['psi1': it[2]] }

            ch_suppa_local_contrasts_psi2 = ch_suppa_local_contrasts_psi1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_local_psi_conditions.for_psi2, by: 0 )
                .map { it -> it[1] + ['psi2': it[2]] }

            // Create input channels to diffsplice process

            ch_split_tpms_events_psi = ch_suppa_local_contrasts_psi2
                .map { it ->
                    [
                        [ id: 'local_' + it.treatment + "-" + it.control ],
                        it.treatment, it.control, it.tpm1, it.tpm2, it.psi1, it.psi2
                    ]
                }

            //
            // Update meta.id of ioe events channel to match the contrasts channel
            //
            ch_events_ioe_for_diffsplice = ch_split_tpms_events_psi
                .map { meta, _treatment, _control, _tpm1, _tpm2, _psi1, _psi2 -> meta }
                .combine(
                    ch_events_ioe
                        .map { _meta, ioe -> ioe }
                )

            DIFFSPLICE_IOE(
                ch_events_ioe_for_diffsplice,
                ch_split_tpms_events_psi.map { meta, treatment, _control, tpm1, _tpm2, psi1, _psi2 -> [ meta, treatment, tpm1, psi1 ] },
                ch_split_tpms_events_psi.map { meta, _treatment, control, _tpm1, tpm2, _psi1, psi2 -> [ meta, control, tpm2, psi2 ] },
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
                diffsplice_nan_threshold
            )

            ch_events_dpsi     = DIFFSPLICE_IOE.out.dpsi
            ch_events_psivec   = DIFFSPLICE_IOE.out.psivec

            if (clusterevents_local_event) {

                // Get ranges for cluster analysis

                ch_events_dpsi_with_conditions = ch_events_dpsi
                    .join(ch_split_tpms_events_psi, by: 0)
                    .map {
                        meta, dpsi, treatment, control, _tpm1, _tpm2, _psi1, _psi2 -> [ meta, treatment, control, dpsi ]
                    }
                ch_events_psivec_with_conditions = ch_events_psivec
                    .join(ch_split_tpms_events_psi, by: 0)
                    .map {
                        meta, psivec, treatment, control, _tpm1, _tpm2, _psi1, _psi2 -> [ meta, treatment, control, psivec ]
                    }

                CLUSTERGROUPS_IOE ( ch_events_psivec_with_conditions )

                ch_groups_ioe = CLUSTERGROUPS_IOE.out.groups

                // Join channels to ensure consistent order

                ch_clusterevents_ioe = ch_events_dpsi_with_conditions
                    .join(ch_events_psivec_with_conditions, by: [0, 1, 2])
                    .join(ch_groups_ioe, by: [0, 1, 2])
                    .map {
                        _meta, cond1, cond2, dpsi, psivec, groups ->
                            [ [ id: 'local_' + cond1 + "-" + cond2 ], dpsi, psivec, groups ]
                    }

                // Run Clustering

                CLUSTEREVENTS_IOE(
                    ch_clusterevents_ioe.map { meta, dpsi, psivec, _groups -> [ meta, dpsi, psivec ] },
                    clusterevents_sigthreshold,
                    clusterevents_dpsithreshold,
                    clusterevents_eps,
                    clusterevents_metric,
                    clusterevents_separation,
                    clusterevents_min_pts,
                    ch_clusterevents_ioe.map { _meta, _dpsi, _psivec, groups -> groups.text.trim() },
                    clusterevents_method
                )

                ch_events_clustvec = CLUSTEREVENTS_IOE.out.clustvec
            }
        }
    }

    // If per isoform analysis:

    ch_events_ioi        = channel.empty()
    ch_isoform_psi       = channel.empty()
    ch_split_isoform_psi = channel.empty()

    ch_isoform_dpsi      = channel.empty()
    ch_isoform_psivec    = channel.empty()

    ch_groups_ioi        = channel.empty()
    ch_isoform_clustvec  = channel.empty()

    if (suppa_per_isoform) {

        // Generate events - transcript level

        SUPPA_GENERATEEVENTS_IOI (
            ch_gtf,
            'ioi',
            generateevents_pool_genes,
            generateevents_event_type,
            generateevents_boundary,
            generateevents_threshold,
            generateevents_exon_length
        )

        ch_events_ioi = SUPPA_GENERATEEVENTS_IOI.out.events

        // Get the psi values per isoform

        SUPPA_PSIPERISOFORM (ch_tpm, ch_gtf)

        ch_isoform_psi = SUPPA_PSIPERISOFORM.out.psi

        // Split the PSI files between the conditions

        SPLIT_FILES_IOI (
            ch_isoform_psi,
            ch_samplesheet,
            '.psi',
            true
        )

        ch_split_isoform_psi = SPLIT_FILES_IOI.out.psis

        // Calculate differential analysis between conditions - Transcript level

        if (diffsplice_transcript_event) {

            // Create contrasts channel

            ch_suppa_isoform_contrasts_raw = ch_contrastsheet

            // Add TPM files to contrasts channel

            SPLIT_FILES_TPM.out.tpms
                .flatten()
                .map { it -> [it.baseName, it ] }
                .multiMap { base, file ->
                    for_tpm1: [ base, file ]
                    for_tpm2: [ base, file ]
                }
                .set { ch_suppa_isoform_tpm_conditions }

            ch_suppa_isoform_contrasts_tpm1 = ch_suppa_isoform_contrasts_raw
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_isoform_tpm_conditions.for_tpm1, by: 0)
                .map { it -> it[1] + ['tpm1': it[2]] }

            ch_suppa_isoform_contrasts_tpm2 = ch_suppa_isoform_contrasts_tpm1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_isoform_tpm_conditions.for_tpm2, by: 0)
                .map { it -> it[1] + ['tpm2': it[2]] }

            // Add PSI files to contrasts channel

            SPLIT_FILES_IOI.out.psis
                .flatten()
                .map { it -> [ it.baseName.toString().replaceAll('transcript_' , ''), it ] }
                .multiMap { base, file ->
                    for_psi1: [ base, file ]
                    for_psi2: [ base, file ]
                }
                .set { ch_suppa_isoform_psi_conditions }

            ch_suppa_isoform_contrasts_psi1 = ch_suppa_isoform_contrasts_tpm2
                .map { it -> [it['treatment'], it] }
                .combine ( ch_suppa_isoform_psi_conditions.for_psi1, by: 0 )
                .map { it -> it[1] + ['psi1': it[2]] }

            ch_suppa_isoform_contrasts_psi2 = ch_suppa_isoform_contrasts_psi1
                .map { it -> [it['control'], it] }
                .combine ( ch_suppa_isoform_psi_conditions.for_psi2, by: 0 )
                .map { it -> it[1] + ['psi2': it[2]] }

            // Create input channels to diffsplice process

            ch_split_tpms_isoform_psi = ch_suppa_isoform_contrasts_psi2
                .map { it ->
                    [
                        [ id: 'transcript_' + it.treatment + "-" + it.control ],
                        it.treatment, it.control, it.tpm1, it.tpm2, it.psi1, it.psi2
                    ]
                }

            //
            // Update meta.id of ioi events channel to match the contrasts channel
            //
            ch_events_ioi_for_diffsplice = ch_split_tpms_isoform_psi
                .map { meta, _treatment, _control, _tpm1, _tpm2, _psi1, _psi2 -> meta }
                .combine(
                    ch_events_ioi
                        .map { _meta, ioi -> ioi }
                )

            DIFFSPLICE_IOI(
                ch_events_ioi_for_diffsplice,
                ch_split_tpms_isoform_psi.map { meta, treatment, _control, tpm1, _tpm2, psi1, _psi2 -> [ meta, treatment, tpm1, psi1 ] },
                ch_split_tpms_isoform_psi.map { meta, _treatment, control, _tpm1, tpm2, _psi1, psi2 -> [ meta, control, tpm2, psi2 ] },
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

            ch_isoform_dpsi   = DIFFSPLICE_IOI.out.dpsi
            ch_isoform_psivec = DIFFSPLICE_IOI.out.psivec

            if (clusterevents_transcript_event) {

                // Get ranges for cluster analysis

                ch_events_dpsi_with_conditions = ch_isoform_dpsi
                    .join(ch_split_tpms_isoform_psi, by: 0)
                    .map {
                        meta, dpsi, treatment, control, _tpm1, _tpm2, _psi1, _psi2 -> [ meta, treatment, control, dpsi ]
                    }

                ch_events_psivec_with_conditions = ch_isoform_psivec
                    .join(ch_split_tpms_isoform_psi, by: 0)
                    .map {
                        meta, psivec, treatment, control, _tpm1, _tpm2, _psi1, _psi2 -> [ meta, treatment, control, psivec ]
                    }

                CLUSTERGROUPS_IOI ( ch_events_psivec_with_conditions )

                ch_groups_ioi = CLUSTERGROUPS_IOI.out.groups

                // Join channels to ensure consistent order

                ch_clusterevents_ioi = ch_events_dpsi_with_conditions
                    .join(ch_events_psivec_with_conditions, by: [0, 1, 2])
                    .join(ch_groups_ioi, by: [0, 1, 2])
                    .map {
                        _meta, cond1, cond2, dpsi, psivec, groups ->
                            [ [ id: 'transcript_' + cond1 + "-" + cond2 ], dpsi, psivec, groups ]
                    }

                // Run Clustering

                CLUSTEREVENTS_IOI(
                    ch_clusterevents_ioi.map { meta, dpsi, psivec, _groups -> [ meta, dpsi, psivec ] },
                    clusterevents_sigthreshold,
                    clusterevents_dpsithreshold,
                    clusterevents_eps,
                    clusterevents_metric,
                    clusterevents_separation,
                    clusterevents_min_pts,
                    ch_clusterevents_ioi.map { _meta, _dpsi, _psivec, groups -> groups.text.trim() },
                    clusterevents_method,
                )

                ch_isoform_clustvec = CLUSTEREVENTS_IOI.out.clustvec
            }
        }
    }

    // Define output

    emit:

    ioe_events              = ch_events_ioe             //    path: events ioe
    ioi_events              = ch_events_ioi             //    path: events ioi

    suppa_local_psi         = ch_events_psi             //    path: suppa_local.psi
    suppa_isoform_psi       = ch_isoform_psi            //    path: suppa_isoform.psi

    split_suppa_tpms        = ch_split_tpms             //    path: suppa_cond1.tpm, suppa_cond2.tpm
    split_suppa_local_psi   = ch_split_events_psi       //    path: suppa_local_cond1.psi, suppa_local_cond2.psi
    split_suppa_isoform_psi = ch_split_isoform_psi      //    path: suppa_isoform_cond1.psi, suppa_isoform_cond2.psi

    dpsi_local              = ch_events_dpsi            //    path: local.dpsi
    psivec_local            = ch_events_psivec          //    path: local.psivec
    dpsi_isoform            = ch_isoform_dpsi           //    path: isoform.dpsi
    psivec_isoform          = ch_isoform_psivec         //    path: isoform.psivec

    cluster_vec_local       = ch_events_clustvec        //    path: local.clustvec
    cluster_vec_isoform     = ch_isoform_clustvec       //    path: isoform.clustvec
}
