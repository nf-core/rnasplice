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

    // If per AS local analysis:

    if (params.suppa_per_local_event) {

        file_type = 'ioe'

        // Generate AS events on the GTF - Local events

        IOE (
            ch_gtf,
            file_type
        )

        // Calculate the psi values of Local events (using events file and TPM)

        PSIPEREVENT (
            IOE.out.events,
            ch_tpm
        )

        // Split the PSI files between the conditions

        output_type = ".psi"
        calc_ranges = true
        prefix = "local"

        SPLIT_PSI_IOE (
            PSIPEREVENT.out.psi,
            ch_samplesheet,
            output_type,
            calc_ranges,
            prefix
        )

        // Calculate differential analysis between conditions

        DIFFSPLICE_IOE(
            IOE.out.events,
            SPLIT_TPM.out.tpms,
            SPLIT_PSI_IOE.out.psis,
            prefix
        )

        // Get ranges for cluster analysis

        SPLIT_PSI_IOE.out.ranges.splitText( by: 1 ){ it.trim() }.set{ ch_ranges_ioe }

        // Run Clustering

        CLUSTEREVENTS_IOE(
            DIFFSPLICE_IOE.out.dpsi,
            DIFFSPLICE_IOE.out.psivec,
            ch_ranges_ioe,
            prefix
        )

    }

    // If per isoform analysis:

    if (params.suppa_per_isoform) {

        file_type = 'ioi'

        // Generate events - transcript level

        IOI (
            ch_gtf,
            file_type
        )

        // Get the psi values per isoform

        PSIPERISOFORM (
            ch_gtf,
            ch_tpm
        )

        // Split the PSI files between the conditions

        output_type = ".psi"
        calc_ranges = true
        prefix = "transcript"

        SPLIT_PSI_IOI (
            PSIPERISOFORM.out.psi,
            ch_samplesheet,
            output_type,
            calc_ranges,
            prefix
        )

        // Calculate differential analysis between conditions - Transcript level

        DIFFSPLICE_IOI(
            IOI.out.events,
            SPLIT_TPM.out.tpms,
            SPLIT_PSI_IOI.out.psis,
            prefix
        )

        // Get ranges for cluster analysis

        SPLIT_PSI_IOI.out.ranges.splitText( by: 1 ){ it.trim() }.set{ ch_ranges_ioi }

        // Run Clustering

        CLUSTEREVENTS_IOI(
            DIFFSPLICE_IOI.out.dpsi,
            DIFFSPLICE_IOI.out.psivec,
            ch_ranges_ioi,
            prefix
        )
    }

    // Define output

    emit:

    ioe_events              = IOE.out.events                     //    path: events ioe
    ioi_events              = IOI.out.events                     //    path: events ioi

    suppa_local_psi         = PSIPEREVENT.out.psi                //    path: suppa_local.psi
    suppa_isoform_psi       = PSIPERISOFORM.out.psi              //    path: suppa_isoform.psi

    split_suppa_tpms        = SPLIT_TPM.out.tpms                 //    path: suppa_cond1.tpm, suppa_cond2.tpm
    split_suppa_local_psi   = SPLIT_PSI_IOE.out.psis             //    path: suppa_local_cond1.psi, suppa_local_cond2.psi
    split_suppa_isoform_psi = SPLIT_PSI_IOI.out.psis             //    path: suppa_isoform_cond1.psi, suppa_isoform_cond2.psi

    dpsi_local              = DIFFSPLICE_IOE.out.dpsi            //    path: local.dpsi
    psivec_local            = DIFFSPLICE_IOE.out.psivec          //    path: local.psivec
    dpsi_isoform            = DIFFSPLICE_IOI.out.dpsi            //    path: isoform.dpsi
    psivec_isoform          = DIFFSPLICE_IOI.out.psivec          //    path: isoform.psivec

    cluster_vec_local       = CLUSTEREVENTS_IOE.out.clustvec     //    path: local.clustvec
    cluster_log_local       = CLUSTEREVENTS_IOE.out.cluster_log  //    path: local.log
    cluster_vec_isoform     = CLUSTEREVENTS_IOI.out.clustvec     //    path: isoform.clustvec
    cluster_log_isoform     = CLUSTEREVENTS_IOI.out.cluster_log  //    path: isoform.log

    versions = ch_versions.ifEmpty(null)                         // channel: [ versions.yml ]
}

