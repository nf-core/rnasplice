# Metro map

The pipeline overview metro map is generated from `assets/metro_map.mmd` using [nf-metro](https://github.com/seqeralabs/nf-metro) 3.0.0 or later. See the [nf-core workflow schematics guide](https://nf-co.re/docs/community/brand/workflow-schematics) for `%%metro` directive syntax and general usage.

If you add or rename pipeline steps, update the `.mmd` source and regenerate the images. The `%%metro output:` directives in the source list every image, so one command writes them all to `docs/images/`:

```bash
pip install 'nf-metro>=3.0.0'
nf-metro render assets/metro_map.mmd
```

This writes:

- `nf-core-rnasplice_metro_map.svg`, which follows the viewer's light or dark mode and is the `manifest.diagram`
- `nf-core-rnasplice_metro_map_animated.svg`, shown in the README
- `nf-core-rnasplice_metro_map_light.png` and `nf-core-rnasplice_metro_map_dark.png`

## How the map is drawn

The lines are the analysis families, so a station sits on a line when its module processes data for that family. Every input type the samplesheet accepts enters at its first consumer: FASTQ files at the start, genome BAM files at SAMtools, transcriptome BAM files at Salmon and Salmon results at tximport. The GTF icons in the QC and exon sections stand for the annotation those tools read. They also give those sections the same top edge as the pre-processing section, whose genome BAM icon sits above its trunk.

Some routes are simplified:

- The Salmon station covers Salmon on the STAR transcriptome BAM files (`--aligner star_salmon`) and Salmon pseudo-alignment of the reads (`--pseudo_aligner salmon`). Both run by default.
- The DTU and event lines reach Salmon through SAMtools, though Salmon reads the STAR transcriptome BAM files and the reads.
- MultiQC only sits on the QC line, though it also reports the Salmon, featureCounts and LeafCutter strand summaries.
- IsoformSwitchAnalyzeR follows tximport, though it reads the Salmon results directly.

The gene annotation tables (`GENEANNOTATION_*`, `GFFREAD_TX2GENE`) and the BED12 annotation for RSeQC (`EAUTILS_GTF2BED`) have no station.

## Live progress overlay

`%%metro process:` directives tie each station to its Nextflow processes, so `nf-metro serve` can light stations up as a run progresses:

```bash
nf-metro serve assets/metro_map.mmd --open --shutdown-after-complete -- \
    nextflow run nf-core/rnasplice -profile test,docker --outdir results
```

After editing the map, check the mapping against the processes of a real run. The mappings are scoped to `NFCORE_RNASPLICE:RNASPLICE`, so they need the full process names of a trace file: `--dag` reads only the short names.

```bash
nextflow run nf-core/rnasplice -profile test,docker --outdir results -with-trace trace.txt
awk -F'\t' 'NR > 1 { sub(/ \(.*/, "", $4); print $4 }' trace.txt | sort -u > processes.txt
nf-metro check-mapping assets/metro_map.mmd --processes processes.txt \
    --ignore 'GENEANNOTATION_' --ignore 'GFFREAD_TX2GENE' --ignore 'EAUTILS_GTF2BED'
```

The test profile does not run every input route. The genome BAM, transcriptome BAM and Salmon results routes run with the `test_genome_bam`, `test_transcriptome_bam` and `test_salmon_results` profiles, and LeafCutter, ASpli, the bigWig tracks and the RSeQC strand check need `--leafcutter --aspli --skip_bigwig false --strand_check`.
