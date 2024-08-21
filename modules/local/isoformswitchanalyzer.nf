process ISOFORMSWITCHANALYZER {
    label 'process_medium'

    conda "bioconda::bioconductor-isoformswitchanalyzer==2.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-isoformswitchanalyzer:2.2.0--r43ha9d7317_0' }"

    input:
    path salmon_output
    path gtf                    // path to gtf file
    path transcript_sequences   // path to isoform nt sequences fasta
    path samplesheet            // path samplesheet
    path contrastsheet          // path contrastsheet
    val  alpha                  // alpha value for differential isoform expression
    val  dIF                    // dIF cutoff value for differential isoform expression

    output:
    path "isoformswitchanalyzer_summary.csv"            , emit: isoformswitchanalyzer_summary
    path "isoformswitchanalyzer_isoformfeatures.csv"    , emit: isoformswitchanalyzer_isoformFeatures
    path "switchlist.rds"                               , emit: switchlist_rds
    path "versions.yml"                                 , emit: versions
    path "results"                                      , emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p results

    run_isoformswitchanalyzer.R $gtf $transcript_sequences $samplesheet $contrastsheet $alpha $dIF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-stager: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
    END_VERSIONS
    """
}
