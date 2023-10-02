//
// Check input contrastsheet and get read channels
//

include { CONTRASTSHEET_CHECK } from '../../modules/local/contrastsheet_check'

workflow CONTRASTS_CHECK {

    take:
    contrastsheet // file: /path/to/contrastsheet.csv

    main:
    CONTRASTSHEET_CHECK ( contrastsheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .set { contrasts }


    emit:
    contrasts
    versions = CONTRASTSHEET_CHECK.out.versions // channel: [ versions.yml ]

}
