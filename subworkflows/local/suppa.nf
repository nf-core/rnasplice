//
// SUPPA Subworkflow
//

include { GENERATE_EVENTS as IOE } from '../../modules/local/suppa_generateevents.nf'
include { GENERATE_EVENTS as IOI } from '../../modules/local/suppa_generateevents.nf'

include { PSIPEREVENT   } from '../../modules/local/suppa_psiperevent.nf'
include { PSIPERISOFORM } from '../../modules/local/suppa_psiperisoform.nf'

include { SPLIT_FILES as SPLIT_TPM  } from '../../modules/local/suppa_split_files.nf'
include { SPLIT_FILES as SPLIT_PSI_IOE  } from '../../modules/local/suppa_split_files.nf'
include { SPLIT_FILES as SPLIT_PSI_IOI  } from '../../modules/local/suppa_split_files.nf'

include { DIFFSPLICE as DIFFSPLICE_IOE } from '../../modules/local/suppa_diffsplice.nf'
include { DIFFSPLICE as DIFFSPLICE_IOI } from '../../modules/local/suppa_diffsplice.nf'

include { CLUSTEREVENTS as CLUSTEREVENTS_IOE } from '../../modules/local/suppa_clusterevents.nf'
include { CLUSTEREVENTS as CLUSTEREVENTS_IOI } from '../../modules/local/suppa_clusterevents.nf'

workflow SUPPA {

    take:

    ch_gtf
    ch_tpm
    ch_samplesheet

    main:

    // define empty versions channel

    ch_versions = Channel.empty()

    // Split the tpm file (contains all samples) into individual files based on condition

    def output_type = ".tpm"
    def calc_ranges = false
    def prefix = ''
    def file_type = ''

    SPLIT_TPM (
        ch_tpm,
        ch_samplesheet,
        output_type,
        calc_ranges,
        prefix
    )

    ch_split_suppa_tpms = SPLIT_TPM.out.tpms

    // If per AS local analysis:

    ch_ioe_events             = Channel.empty()
    ch_suppa_local_psi        = Channel.empty()
    ch_split_suppa_local_psi  = Channel.empty()

    ch_dpsi_local             = Channel.empty()
    ch_psivec_local           = Channel.empty()

    ch_cluster_vec_local      = Channel.empty()
    ch_cluster_log_local      = Channel.empty()

    if (params.suppa_per_local_event) {

        file_type = 'ioe'

        // Generate AS events on the GTF - Local events

        IOE (
            ch_gtf,
            file_type
        )

        ch_ioe_events  = IOE.out.events

        // Calculate the psi values of Local events (using events file and TPM)

        PSIPEREVENT (
            ch_ioe_events,
            ch_tpm
        )

        ch_suppa_local_psi = PSIPEREVENT.out.psi

        // Split the PSI files between the conditions

        output_type = ".psi"
        calc_ranges = true
        prefix = "local"

        SPLIT_PSI_IOE (
            ch_suppa_local_psi,
            ch_samplesheet,
            output_type,
            calc_ranges,
            prefix
        )

        ch_split_suppa_local_psi = SPLIT_PSI_IOE.out.psis

        // Calculate differential analysis between conditions

        if (params.diffsplice_local_event) {

            DIFFSPLICE_IOE(
                ch_ioe_events,
                ch_split_suppa_tpms,
                ch_split_suppa_local_psi,
                prefix
            )

            ch_dpsi_local     = DIFFSPLICE_IOE.out.dpsi
            ch_psivec_local   = DIFFSPLICE_IOE.out.psivec

            if (params.clusterevents_local_event) {

                // Get ranges for cluster analysis

                SPLIT_PSI_IOE.out.ranges.splitText( by: 1 ){ it.trim() }.set{ ch_ranges_ioe }

                // Run Clustering

                CLUSTEREVENTS_IOE(
                    ch_dpsi_local,
                    ch_psivec_local,
                    ch_ranges_ioe,
                    prefix
                )

                ch_cluster_vec_local   = CLUSTEREVENTS_IOE.out.clustvec
                ch_cluster_log_local   = CLUSTEREVENTS_IOE.out.cluster_log
            }
        }
    }

    // If per isoform analysis:

    ch_ioi_events              = Channel.empty()
    ch_suppa_isoform_psi       = Channel.empty()
    ch_split_suppa_isoform_psi = Channel.empty()

    ch_dpsi_isoform            = Channel.empty()
    ch_psivec_isoform          = Channel.empty()

    ch_cluster_vec_isoform     = Channel.empty()
    ch_cluster_log_isoform     = Channel.empty()

    if (params.suppa_per_isoform) {

        file_type = 'ioi'

        // Generate events - transcript level

        IOI (
            ch_gtf,
            file_type
        )

        ch_ioi_events = IOI.out.events

        // Get the psi values per isoform

        PSIPERISOFORM (
            ch_gtf,
            ch_tpm
        )

        ch_suppa_isoform_psi = PSIPERISOFORM.out.psi

        // Split the PSI files between the conditions

        output_type = ".psi"
        calc_ranges = true
        prefix = "transcript"

        SPLIT_PSI_IOI (
            ch_suppa_isoform_psi,
            ch_samplesheet,
            output_type,
            calc_ranges,
            prefix
        )

        ch_split_suppa_isoform_psi = SPLIT_PSI_IOI.out.psis

        // Calculate differential analysis between conditions - Transcript level

        if (params.diffsplice_isoform) {

            DIFFSPLICE_IOI(
                ch_ioi_events,
                ch_split_suppa_tpms,
                ch_split_suppa_isoform_psi,
                prefix
            )

            ch_dpsi_isoform   = DIFFSPLICE_IOI.out.dpsi
            ch_psivec_isoform = DIFFSPLICE_IOI.out.psivec

            if (params.clusterevents_isoform) {

                // Get ranges for cluster analysis

                SPLIT_PSI_IOI.out.ranges.splitText( by: 1 ){ it.trim() }.set{ ch_ranges_ioi }

                // Run Clustering

                CLUSTEREVENTS_IOI(
                    ch_dpsi_isoform,
                    ch_psivec_isoform,
                    ch_ranges_ioi,
                    prefix
                )

                ch_cluster_vec_isoform = CLUSTEREVENTS_IOI.out.clustvec
                ch_cluster_log_isoform = CLUSTEREVENTS_IOI.out.cluster_log
            }
        }
    }

    // Define output

    emit:

    ioe_events              = ch_ioe_events                      //    path: events ioe
    ioi_events              = ch_ioi_events                      //    path: events ioi

    suppa_local_psi         = ch_suppa_local_psi                 //    path: suppa_local.psi
    suppa_isoform_psi       = ch_suppa_isoform_psi               //    path: suppa_isoform.psi

    split_suppa_tpms        = ch_split_suppa_tpms                //    path: suppa_cond1.tpm, suppa_cond2.tpm
    split_suppa_local_psi   = ch_split_suppa_local_psi           //    path: suppa_local_cond1.psi, suppa_local_cond2.psi
    split_suppa_isoform_psi = ch_split_suppa_isoform_psi         //    path: suppa_isoform_cond1.psi, suppa_isoform_cond2.psi

    dpsi_local              = ch_dpsi_local                      //    path: local.dpsi
    psivec_local            = ch_psivec_local                    //    path: local.psivec
    dpsi_isoform            = ch_dpsi_isoform                    //    path: isoform.dpsi
    psivec_isoform          = ch_psivec_isoform                  //    path: isoform.psivec

    cluster_vec_local       = ch_cluster_vec_local               //    path: local.clustvec
    cluster_log_local       = ch_cluster_log_local               //    path: local.log
    cluster_vec_isoform     = ch_cluster_vec_isoform             //    path: isoform.clustvec
    cluster_log_isoform     = ch_cluster_log_isoform             //    path: isoform.log

    versions = ch_versions.ifEmpty(null)                         // channel: [ versions.yml ]
}

