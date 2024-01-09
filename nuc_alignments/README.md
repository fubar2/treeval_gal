[nuc_alignments](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/nuc_alignments.nf) is the nucleotide pathway for hic_mapping

This subworkflow is used in [hic_mapping](hic_mapping) as the nucleotide aligner.
![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_gene_alignment.png)

This DDL has function calls explained below.
Most of the rest of the DDL is not going to be needed other than to
figure out exactly how each function gets parameters supplied to the actual command lines.
Tools needed to handle the PAF format can easily be avoided - the Galaxy minimap2 tool can create a bam output automatically.

```

workflow NUC_ALIGNMENTS {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    reference_index     // Channel: tuple [ val(meta), path(file) ]
    nuc_files           // Channel: tuple [ val(meta), path(file) ]
    dot_genome          // Channel: tuple [ val(meta), path(file) ]
    intron_size         // Channel: val(50k)

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: COLLECTION FROM GENE_ALIGNMENT IS A LIST OF ALL META AND ALL FILES
    //        BELOW CONVERTS INTO TUPLE FORMAT AND ADDS BOOLEANS FOR MINIMAP2_ALIGN
    //
    nuc_files
        .flatten()
        .buffer( size: 2 )
        .combine ( reference_tuple )
        .combine( intron_size )
        .map { meta, nuc_file, ref_meta, ref, intron ->
            tuple( [id:             meta.id,
                    type:           meta.type,
                    org:            meta.org,
                    intron_size:    intron,
                    split_prefix:   nuc_file.toString().split('/')[-1].split('.fasta')[0],
                    single_end:     true
                    ],
                    nuc_file,
                    ref,
                    true,
                    false,
                    false
            )
        }
        .multiMap { meta, nuc_file, reference, bool_1, bool_2, bool_3 ->
            nuc             : tuple( meta, nuc_file)
            ref             : reference
            bool_bam_output : bool_1
            bool_cigar_paf  : bool_2
            bool_cigar_bam  : bool_3
        }
        .set { formatted_input }

    //
    // MODULE: ALIGNS REFERENCE FAIDX TO THE GENE_ALIGNMENT QUERY FILE FROM NUC_FILES
    //         EMITS ALIGNED BAM FILE
    //
    MINIMAP2_ALIGN (
        formatted_input.nuc,
        formatted_input.ref,
        formatted_input.bool_bam_output,
        formatted_input.bool_cigar_paf,
        formatted_input.bool_cigar_bam
    )
    ch_versions     = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // LOGIC: CONVERTS THE MINIMAP OUTPUT TUPLE INTO A GROUPED TUPLE PER INPUT QUERY ORGANISM
    //        AND DATA TYPE (RNA, CDS, DNA).
    //
    MINIMAP2_ALIGN.out.bam
        .map { meta, file ->
            tuple(
                [   id: meta.org,
                    type: meta.type ],
                file) }
        .groupTuple( by: [0] )  // group by meta list
        .set { merge_input }

    //
    // MODULE: MERGES THE BAM FILES FOUND IN THE GROUPED TUPLE IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM
    SAMTOOLS_MERGE (
        merge_input,
        reference_tuple,
        reference_index
    )
    ch_versions     = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // SUBWORKFLOW: GENERATES A PUNCHLIST FROM MERGED BAM FILE
    //
    PUNCHLIST (
        reference_tuple,
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions     = ch_versions.mix(PUNCHLIST.out.versions)

    //
    // MODULE: CONVERTS THE ABOVE MERGED BAM INTO BED FORMAT
    //
    BEDTOOLS_BAMTOBED ( SAMTOOLS_MERGE.out.bam )
    ch_versions     = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    // TODO: try filtering out here too

    //
    // MODULE: SORTS THE ABOVE BED FILE
    //
    BEDTOOLS_SORT (
        BEDTOOLS_BAMTOBED.out.bed,
        []
    )
    ch_versions     = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    //
    // LOGIC: COMBINES GENOME_FILE CHANNEL AND ABOVE OUTPUT, SPLITS INTO TWO CHANNELS
    //        ALSO FILTERS OUT EMPTY MERGED.BED BASED ON WHETHER FILE IS >141 BYTES
    //
    BEDTOOLS_SORT.out.sorted
        .map { meta, file ->
                tuple( [    id:         meta.id,
                            type:       meta.type,
                            file_size:  file.size()
                        ],
                        file ) }
        .filter { it[0].file_size >= 141 } // Take the first item in input (meta) and check if size is more than a symlink
        .combine( dot_genome )
        .multiMap { meta, ref, genome_meta, genome ->
            bed_file:   tuple( [    id:         meta.id,
                                    type:       meta.type,
                                ],
                                ref )
            dot_genome: genome
        }
        .set { ucsc_input }

    //
    // MODULE: CONVERTS GENOME FILE AND BED INTO A BIGBED FILE
    //
    UCSC_BEDTOBIGBED (
        ucsc_input.bed_file,
        ucsc_input.dot_genome,
        []
    )
```

[MINIMAP2_ALIGN](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/minimap2/align/main.nf)
Uses *conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"* to run

```
 minimap2 \\
        $args \\
        -t $task.cpus \\
        "${reference ?: reads}" \\
        "$reads" \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output
```

[SAMTOOLS_MERGE](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/samtools/merge/main.nf) uses *conda "bioconda::samtools=1.17"* to run

```
samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${reference} \\
        ${prefix}.${file_type} \\
        $input_files
```

[PUNCHLIST](../punchlist) is another subworkflow

[BEDTOOLS_BAMTOBED](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/bamtobed/main.nf) calls

```
 bedtools \\
        bamtobed \\
        $args \\
        -i $bam \\
        > ${prefix}.bed
```

bedtools sort is then needed and available
followed by [UCSC-bedtobigbed](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/ucsc-bedtobigbed) using:

```
  bedToBigBed
        $bed \\
        $sizes \\
        $as_option \\
        $args \\
        ${prefix}.bigBed
```
Also used in [ancestral_gene](ancestral_gene)


