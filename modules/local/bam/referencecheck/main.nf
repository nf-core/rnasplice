process BAM_REFERENCECHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data':
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path(bam, includeInputs: true), emit: bam
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -H ${bam} \\
        | awk -F '\\t' '\$1 == "@SQ" { name = ""; length_ = ""; for (i = 2; i <= NF; i++) { if (\$i ~ /^SN:/) name = substr(\$i, 4); if (\$i ~ /^LN:/) length_ = substr(\$i, 4) } print name "\\t" length_ }' \\
        > bam_sequences.tsv

    awk -F '\\t' 'NR == FNR { reference[\$1] = \$2; next } !(\$1 in reference) || reference[\$1] != \$2' ${fai} bam_sequences.tsv \\
        > mismatched_sequences.tsv

    if [ ! -s bam_sequences.tsv ]; then
        echo "ERROR: ${bam} has no @SQ lines in its header, so it cannot be checked against ${fai.baseName}." >&2
        exit 1
    fi

    if [ -s mismatched_sequences.tsv ]; then
        {
            echo "ERROR: ${bam} was not aligned to ${fai.baseName}."
            echo "\$(wc -l < mismatched_sequences.tsv) of its \$(wc -l < bam_sequences.tsv) reference sequences are missing from it, or have another length, for example (name, length):"
            head -n 10 mismatched_sequences.tsv
            echo "Give the pipeline the references the BAM files were aligned to. For BAM files from nf-core/rnaseq that is its genome FASTA and, for transcriptome BAM files, the transcript FASTA it wrote to '<outdir>/genome', given with '--transcript_fasta'."
        } >&2
        exit 1
    fi
    """

    stub:
    """
    """
}
