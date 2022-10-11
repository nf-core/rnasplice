//
// SUPPA
//
include { GENERATEEVENTS as IOE  } from '../../modules/local/suppa_generateevents.nf'
include { GENERATEEVENTS  as IOI } from '../../modules/local/suppa_generateevents.nf'

include { PSIPEREVENT  as IOE_PSI } from '../../modules/local/suppa_psiperevent.nf'
include { PSIPEREVENT  as IOI_PSI } from '../../modules/local/suppa_psiperevent.nf'

include { PSIPERISOFORM  as ISOFORM_PSI } from '../../modules/local/suppa_psiperisoform.nf'

include { SPLIT_FILES as SPLIT_TPM  } from '../../modules/local/suppa_split_files.nf'
include { SPLIT_FILES as SPLIT_PSI_IOE  } from '../../modules/local/suppa_split_files.nf'
include { SPLIT_FILES as SPLIT_PSI_IOI  } from '../../modules/local/suppa_split_files.nf'
//include { SPLIT_FILES as SPLIT_PSI_ISO  } from '../../modules/local/split_files.nf' //SUPPA Check!

include { GETLIST } from '../../modules/local/suppa_getlist.nf'

include { DIFFSPLICE as DIFFSPLICE_IOE } from '../../modules/local/suppa_diffsplice.nf'
include { DIFFSPLICE as DIFFSPLICE_IOI } from '../../modules/local/suppa_diffsplice.nf'
//include { DIFFSPLICE_ISO } from '../../modules/local/diffsplice_iso.nf' //SUPPA Check!

include { CLUSTEREVENTS as CLUSTEREVENTS_IOE } from '../../modules/local/suppa_clusterevents.nf'
include { CLUSTEREVENTS as CLUSTEREVENTS_IOI } from '../../modules/local/suppa_clusterevents.nf'

workflow SUPPA {

    take:

    ch_gtf           
    ch_tpm 
    ch_input        

    main:

    // define empty versions channel
    ch_versions = Channel.empty()
    
    // Calculate the AS events on the GTF
    IOE ( ch_gtf, 'ioe') //Local events
    IOI ( ch_gtf, 'ioi') //Transcript events

    // Calculate the psi values of events
    IOE_PSI (IOE.out.events, ch_tpm, 'ioe') // Get the psi values of the Local events (using events file and TPM)
    IOI_PSI (IOI.out.events, ch_tpm, 'ioi') // Get the psi values of the Transcript events (using events file and TPM)
    
    ISOFORM_PSI (ch_gtf, ch_tpm) // Get psi per transcript Isoform (using GTF and TPM) 
    
    // Get Strings as needed for SPLIT_FILES module
    GETLIST(ch_input) // Get the list of samples in each condition in a text file as per the format needed for SPLIT_FILES module
    tpm_list = GETLIST.out.tpm_list.getText() // Fetch the string in the text file 
    psi_list = GETLIST.out.psi_list.getText() // Fetch the string in the text file 
  
    // Split the PSI and TPM files between the conditions
    SPLIT_TPM (ch_tpm,tpm_list) // Split the tpm file (contains all samples) into individual files based on condition
    SPLIT_PSI_IOE (IOE_PSI.out.psi,psi_list) // Split the Local events psi file into individual files based on condition
    SPLIT_PSI_IOI (IOI_PSI.out.psi,psi_list) // Split the Transcript events psi file into individual files based on condition
    //SPLIT_PSI_ISO (ISOFORM_PSI.out.psi,psi_list) //SUPPA Check!

    // Calculate differential analysis between conditions
    DIFFSPLICE_IOE(IOE.out.events,SPLIT_TPM.out.tpms,SPLIT_PSI_IOE.out.psis) // Local events
    DIFFSPLICE_IOI(IOI.out.events,SPLIT_TPM.out.tpms,SPLIT_PSI_IOI.out.psis) // Transcript events
    //DIFFSPLICE_ISO(SPLIT_TPM.out.tpms,SPLIT_PSI_ISO.out.psis)

    CLUSTEREVENTS_IOE(DIFFSPLICE_IOE.out.dpsi, DIFFSPLICE_IOE.out.psivec)
    //CLUSTEREVENTS_IOI(DIFFSPLICE_IOI.out.dpsi, DIFFSPLICE_IOI.out.psivec) 

    //Define output
    emit:
    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}