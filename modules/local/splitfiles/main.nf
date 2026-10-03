process SPLIT_FILES {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7aa202fe46bc3ec5eb02d21c659ee8d93b1995724712fface3837047691651a2/data':
        'community.wave.seqera.io/library/r-base:4.6.1--e4f1a108384f0df8' }"

    input:
    tuple val(meta), path(tpm_psi), val(samples)
    val extension // either 'tpm' or 'psi'

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: split
    tuple val("${task.process}"), val('r-base'), eval('R --version 2>&1 | head -n 1 | sed "s/^.*version //; s/ .*$//"'), topic: versions, emit: versions_r

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!['tpm', 'psi'].contains(extension)) {
        error("SPLIT_FILES: extension must be 'tpm' or 'psi', got '${extension}'")
    }
    // The samples as an R character vector, in the order of the columns to write
    def samples_r = samples.collect { sample -> "\"${sample}\"" }.join(', ')
    """
    Rscript - <<'EOF'
    samples <- c(${samples_r})

    # The header of a SUPPA matrix has one field less than its rows, so the first
    # column becomes the row names. The values are read as text, so that they are
    # written back as SUPPA wrote them, its 'nan' included
    input_data <- read.csv(
        "${tpm_psi}",
        sep          = "\\t",
        header       = TRUE,
        check.names  = FALSE,
        colClasses   = "character",
        na.strings   = character(0)
    )

    missing <- setdiff(samples, colnames(input_data))
    if (length(missing) > 0) {
        stop("samples missing from ${tpm_psi}: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    write.table(
        input_data[, samples, drop = FALSE],
        file  = "${prefix}.${extension}",
        quote = FALSE,
        sep   = "\\t"
    )
    EOF
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.${extension}
    """
}
