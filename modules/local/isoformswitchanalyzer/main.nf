process ISOFORMSWITCHANALYZER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.12.0--r45hd2fad28_0'
        : 'biocontainers/bioconductor-isoformswitchanalyzer:2.12.0--r45hd2fad28_0'}"

    input:
    path salmon_output          // path: one Salmon quant directory per sample
    path gtf                    // path: /path/to/genes.gtf
    path transcript_sequences   // path: /path/to/transcripts.fa (isoform nucleotide sequences)
    path samplesheet            // path: /path/to/samplesheet.csv
    path contrastsheet          // path: /path/to/contrastsheet.csv
    val  alpha                  // val:  FDR cutoff for the isoform switch test
    val  dIF                    // val:  minimum absolute difference in isoform fraction

    output:
    path "isoformswitchanalyzer_summary.csv"        , emit: isoformswitchanalyzer_summary
    path "isoformswitchanalyzer_isoformfeatures.csv", emit: isoformswitchanalyzer_isoformFeatures
    path "switchlist.rds"                           , emit: switchlist_rds
    path "results"                                  , emit: results
    path "versions.yml"                             , topic: versions, emit: versions_isoformswitchanalyzer

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'run_isoformswitchanalyzer.R'

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}

    mkdir -p results
    touch isoformswitchanalyzer_summary.csv
    touch isoformswitchanalyzer_isoformfeatures.csv
    touch switchlist.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | head -n 1 | sed 's/^.*version //; s/ .*\$//')
        bioconductor-isoformswitchanalyzer: \$(Rscript -e "library(IsoformSwitchAnalyzeR); cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))")
    END_VERSIONS
    """
}
