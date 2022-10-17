//
// SUPPA Subworkflow
//

include { GENERATE_EVENTS as IOE } from '../../modules/local/suppa_generateevents.nf'
include { GENERATE_EVENTS as IOI } from '../../modules/local/suppa_generateevents.nf'

include { PSIPEREVENT as IOE_PSI } from '../../modules/local/suppa_psiperevent.nf'
include { PSIPEREVENT as IOI_PSI } from '../../modules/local/suppa_psiperevent.nf'

include { PSIPERISOFORM as ISOFORM_PSI } from '../../modules/local/suppa_psiperisoform.nf'

include { SPLIT_FILES as SPLIT_TPM  } from '../../modules/local/suppa_split_files.nf'
include { SPLIT_FILES as SPLIT_PSI_IOE  } from '../../modules/local/suppa_split_files.nf'
include { SPLIT_FILES as SPLIT_PSI_IOI  } from '../../modules/local/suppa_split_files.nf'
// include { SPLIT_FILES as SPLIT_PSI_ISO  } from '../../modules/local/split_files.nf' //SUPPA Check!

//include { GETLIST } from '../../modules/local/suppa_getlist.nf'

include { DIFFSPLICE as DIFFSPLICE_IOE } from '../../modules/local/suppa_diffsplice.nf'
include { DIFFSPLICE as DIFFSPLICE_IOI } from '../../modules/local/suppa_diffsplice.nf'
// include { DIFFSPLICE_ISO } from '../../modules/local/diffsplice_iso.nf' //SUPPA Check!

include { CLUSTEREVENTS as CLUSTEREVENTS_IOE } from '../../modules/local/suppa_clusterevents.nf'
include { CLUSTEREVENTS as CLUSTEREVENTS_IOI } from '../../modules/local/suppa_clusterevents.nf'

workflow SUPPA {

    take:

    ch_gtf           
    ch_tpm 
    //ch_input  
    ch_getlist_suppa_tpm  
    ch_getlist_suppa_psi 

    main:

    // Get Strings as needed for SPLIT_FILES module

    //GETLIST( ch_input ) // Get the list of samples in each condition in a text file as per the format needed for SPLIT_FILES module

    //tpm_list = GETLIST.out.tpm_list.getText() // Fetch the string in the text file 
    //psi_list = GETLIST.out.psi_list.getText() // Fetch the string in the text file 
    
    SPLIT_TPM ( ch_tpm, ch_getlist_suppa_tpm ) // Split the tpm file (contains all samples) into individual files based on condition

    // define empty versions channel
    ch_versions = Channel.empty()

    // If per AS local analysis: 

    if (params.suppa_local_as_ioe) {

        def file_type = 'ioe'

        // Calculate the AS events on the GTF - Local events
        IOE ( 
            ch_gtf, 
            file_type 
        )

        // Calculate the psi values of Local events (using events file and TPM)
        IOE_PSI ( 
            IOE.out.events, 
            ch_tpm, 
            file_type 
        )

        // Split the PSI and TPM files between the conditions
        SPLIT_PSI_IOE ( 
            IOE_PSI.out.psi, 
            ch_getlist_suppa_psi 
        )

        // Calculate differential analysis between conditions
        DIFFSPLICE_IOE(
            IOE.out.events,
            SPLIT_TPM.out.tpms,
            SPLIT_PSI_IOE.out.psis
        )
      /*  DIFFSPLICE_IOE
            .out
            .psivec
            .withReader { line = it.readLine()}
            .map{it -> [it.toString()]}
            .set { ch_test }
            ch_test.view()*/

        CLUSTEREVENTS_IOE(
            DIFFSPLICE_IOE.out.dpsi,
            DIFFSPLICE_IOE.out.psivec
        )

    }

    // If per transcript local analysis: 

    if (params.suppa_local_tx_ioi) {
        
        def file_type = 'ioi'

        // Calculate the AS events on the GTF - Transcript events
        IOI ( 
            ch_gtf, 
            file_type 
        )

        // Get the psi values of the Transcript events (using events file and TPM)
        IOI_PSI ( 
            IOI.out.events, 
            ch_tpm, 
            file_type 
        )

        // Split the PSI and TPM files between the conditions
        SPLIT_PSI_IOI ( 
            IOI_PSI.out.psi, 
            ch_getlist_suppa_psi 
        )

        // Calculate differential analysis between conditions - Transcript events
        DIFFSPLICE_IOI(
            IOI.out.events, 
            SPLIT_TPM.out.tpms,
            SPLIT_PSI_IOI.out.psis
        ) 

        CLUSTEREVENTS_IOI(
            DIFFSPLICE_IOI.out.dpsi,
            DIFFSPLICE_IOI.out.psivec
        ) 

    }
    
    // If isoform analysis required

    if (params.suppa_per_isoform) {

        ISOFORM_PSI (ch_gtf, ch_tpm) // Get psi per transcript Isoform (using GTF and TPM) 
        // SPLIT_PSI_ISO (ISOFORM_PSI.out.psi,psi_list) //SUPPA Check!
        // DIFFSPLICE_ISO(SPLIT_TPM.out.tpms,SPLIT_PSI_ISO.out.psis)
    }

    // Define output

    emit:

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}