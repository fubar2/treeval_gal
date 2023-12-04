### [#3 hic_mapping](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/hic_mapping.nf)

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_hic_mapping.png)
(from [https://github.com/sanger-tol/treeval/blob/dev/docs/output.md#busco-analysis](https://github.com/sanger-tol/treeval/blob/dev/docs/output.md#hic_mapping))

```
Output files
    hic_files/
        *_pretext_hr.pretext: High resolution pretext map.
        *_pretext_normal.pretext: Standard resolution pretext map.
        *.mcool: HiC map required for HiGlass

```
This is a lot of steps. Some local code - could need new tools. Juicer? Cooler? cram_filter?

This DDL has function calls explained below. Most of the rest of the DDL is not going to be needed other than to
figure out exactly how each function gets parameters supplied to the actual command lines.

```
workflow HIC_MAPPING {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path( file )      ]
    reference_index     // Channel: tuple [ val(meta), path( file )      ]
    dot_genome          // Channel: tuple [ val(meta), path( datafile )  ]
    hic_reads_path      // Channel: tuple [ val(meta), path( directory ) ]
    assembly_id         // Channel: val( id )
    gap_file            // Channel: tuple [ val(meta), path( file )      ]
    coverage_file       // Channel: tuple [ val(meta), path( file )      ]
    logcoverage_file    // Channel: tuple [ val(meta), path( file )      ]
    telo_file           // Channel: tuple [ val(meta), path( file )      ]
    repeat_density_file // Channel: tuple [ val(meta), path( file )      ]
    workflow_setting    // Channel: val( { RAPID | FULL } )

    main:
    ch_versions         = Channel.empty()

    // COMMENT: 1000bp BIN SIZE INTERVALS FOR CLOAD
    ch_cool_bin         = Channel.of( 1000 )

    //
    // MODULE: Indexing on reference output the folder of indexing files
    //
    BWAMEM2_INDEX (
        reference_tuple
    )
    ch_versions         = ch_versions.mix( BWAMEM2_INDEX.out.versions )

    //
    // LOGIC: make channel of hic reads as input for GENERATE_CRAM_CSV
    //
    reference_tuple
        .combine( hic_reads_path )
        .map { meta, ref, hic_meta, hic_reads_path ->
                tuple(
                    [ id: meta.id, single_end: true],
                    hic_reads_path
                )
        }
        .set { get_reads_input }

    //
    // MODULE: generate a cram csv file containing the required parametres for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV (
        get_reads_input
    )
    ch_versions         = ch_versions.mix( GENERATE_CRAM_CSV.out.versions )

    //
    // LOGIC: organise all parametres into a channel for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV.out.csv
        .splitCsv()
        .combine ( reference_tuple )
        .combine ( BWAMEM2_INDEX.out.index )
        .map{ cram_id, cram_info, ref_id, ref_dir, bwa_id, bwa_path ->
                tuple([
                        id: cram_id.id
                        ],
                    file(cram_info[0]),
                    cram_info[1],
                    cram_info[2],
                    cram_info[3],
                    cram_info[4],
                    cram_info[5],
                    cram_info[6],
                    bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1]
                )
        }
        .set { ch_filtering_input }

    //
    // MODULE: parallel proccessing bwa-mem2 alignment by given interval of containers from cram files
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT (
        ch_filtering_input
    )
    ch_versions         = ch_versions.mix( CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions )

    //
    // LOGIC: PREPARING BAMS FOR MERGE
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam
        .map{ meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [
                id: file[0].toString().split('/')[-1].split('_')[0] + '_' + file[0].toString().split('/')[-1].split('_')[1]
                ],
                file
            )
        }
        .set { collected_files_for_merge }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MERGE (
        collected_files_for_merge,
        reference_tuple,
        reference_index
    )
    ch_versions         = ch_versions.mix ( SAMTOOLS_MERGE.out.versions.first() )

    //
    // LOGIC: PREPARING PRETEXT MAP INPUT
    //
    SAMTOOLS_MERGE.out.bam
        .combine( reference_tuple )
        .multiMap { bam_meta, bam, ref_meta, ref_fa ->
            input_bam:  tuple( [    id: bam_meta.id,
                                    sz: file( bam ).size() ],
                                bam
                        )
            reference:  ref_fa
        }
        .set { pretext_input }

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR LOW RES
    //
    PRETEXTMAP_STANDRD (
        pretext_input.input_bam,
        pretext_input.reference
    )
    ch_versions         = ch_versions.mix( PRETEXTMAP_STANDRD.out.versions )

    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //
    PRETEXT_INGEST_SNDRD (
        PRETEXTMAP_STANDRD.out.pretext,
        gap_file,
        coverage_file,
        logcoverage_file,
        telo_file,
        repeat_density_file
    )
    ch_versions         = ch_versions.mix( PRETEXT_INGEST_SNDRD.out.versions )

    //
    // LOGIC: HIRES IS TOO INTENSIVE FOR RUNNING IN GITHUB CI SO THIS STOPS IT RUNNING
    //
    if ( params.config_profile_name ) {
        config_profile_name = params.config_profile_name
    } else {
        config_profile_name = 'Local'
    }

    if ( !config_profile_name.contains('GitHub') ) {
        //
        // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR HIGH RES
        //
        PRETEXTMAP_HIGHRES (
            pretext_input.input_bam,
            pretext_input.reference
        )
        ch_versions         = ch_versions.mix( PRETEXTMAP_HIGHRES.out.versions )

        //
        // NOTICE: This could fail on LARGE hires maps due to some memory parameter in the C code
        //         of pretext graph. There is a "fixed" version in sanger /software which may need
        //         to be released in this case
        //
        PRETEXT_INGEST_HIRES (
            PRETEXTMAP_HIGHRES.out.pretext,
            gap_file,
            coverage_file,
            logcoverage_file,
            telo_file,
            repeat_density_file
        )
        ch_versions         = ch_versions.mix( PRETEXT_INGEST_HIRES.out.versions )
    }

    //
    // MODULE: GENERATE PNG FROM STANDARD PRETEXT
    //
    SNAPSHOT_SRES (
        PRETEXTMAP_STANDRD.out.pretext
    )
    ch_versions         = ch_versions.mix ( SNAPSHOT_SRES.out.versions )

    // NOTE: CURRENTLY UNDER INVESTIGATION
    //
    // MODULE: GENERATE PNG FROM HIGHRES PRETEXT
    //
    // SNAPSHOT_HRES ( PRETEXTMAP_HIGHRES.out.pretext )
    // ch_versions         = ch_versions.mix ( SNAPSHOT_HRES.out.versions )

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MARKDUP (
        pretext_input.input_bam,
        pretext_input.reference
    )
    ch_versions         = ch_versions.mix ( SAMTOOLS_MARKDUP.out.versions )

    //
    // MODULE: SAMTOOLS FILTER OUT DUPLICATE READS | BAMTOBED | SORT BED FILE
    //
    BAMTOBED_SORT(
        SAMTOOLS_MARKDUP.out.bam
    )
    ch_versions         = ch_versions.mix( BAMTOBED_SORT.out.versions )

    //
    // MODULE: GENERATE CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED( BAMTOBED_SORT.out.sorted_bed )
    ch_versions         = ch_versions.mix( GET_PAIRED_CONTACT_BED.out.versions )

    //
    // LOGIC: SECTION ONLY NEEDED FOR TREEVAL VISUALISATION, NOT RAPID ANALYSIS
    //
    if (workflow_setting == 'FULL' && !config_profile_name.contains('GitHub')) {
        //
        // LOGIC: PREPARE JUICER TOOLS INPUT
        //
        GET_PAIRED_CONTACT_BED.out.bed
            .combine( dot_genome )
            .multiMap {  meta, paired_contacts, meta_my_genome, my_genome ->
                paired      :   tuple([ id: meta.id, single_end: true], paired_contacts )
                genome      :   my_genome
                id          :   meta.id
            }
            .set { ch_juicer_input }

        //
        // MODULE: GENERATE HIC MAP, ONLY IS PIPELINE IS RUNNING ON ENTRY FULL
        //

        JUICER_TOOLS_PRE(
            ch_juicer_input.paired,
            ch_juicer_input.genome,
            ch_juicer_input.id
        )
        ch_versions         = ch_versions.mix( JUICER_TOOLS_PRE.out.versions )
    }

    //
    // LOGIC: BIN CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED.out.bed
        .join( BAMTOBED_SORT.out.sorted_bed )
        .combine( ch_cool_bin )
        .set { ch_binned_pairs }

    //
    // LOGIC: PREPARE COOLER INPUT
    //
    ch_binned_pairs
        .combine(dot_genome)
        .multiMap { meta, pairs, bed, cool_bin, meta_my_genome, my_genome ->
            cooler_in   : tuple ( meta, pairs, bed, cool_bin )
            genome_file : my_genome
        }
        .set { ch_cooler }

    //
    // MODULE: GENERATE A MULTI-RESOLUTION COOLER FILE BY COARSENING
    //
    COOLER_CLOAD(
        ch_cooler.cooler_in,
        ch_cooler.genome_file
    )
    ch_versions         = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // LOGIC: REFACTOR CHANNEL FOR ZOOMIFY
    //
    COOLER_CLOAD.out.cool
        .map{ meta, cools, cool_bin ->
            [meta, cools]
        }
        .set{ch_cool}

    //
    // MODULE: ZOOM COOL TO MCOOL
    //
    COOLER_ZOOMIFY(ch_cool)
    ch_versions         = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    //
    // LOGIC: FOR REPORTING
    //

    ch_cram_files = GrabFiles( get_reads_input )

    ch_cram_files
        .collect()
        .map { meta, cram ->
            tuple( [    id: 'cram',
                        sz: cram instanceof ArrayList ? cram.collect { it.size()} : cram.size(),
                    ],
                    cram
            )
        }
        .combine( GENERATE_CRAM_CSV.out.csv )
        .map { meta, data, meta2, csv ->
            tuple( [    id: meta.id,
                        sz: meta.sz,
                        cn: csv.countLines()
                    ],
                    data
            )
        }
        .set { ch_reporting_cram }
```

## Steps broken down

BWAMEM2_INDEX is an NF-core module that executes:

```
mkdir bwamem2
    bwa-mem2 \\
        index \\
        $args \\
        $fasta -p bwamem2/${fasta}
```

[GENERATE_CRAM_CSV](https://github.com/sanger-tol/treeval/blob/dev/modules/local/generate_cram_csv.nf) uses a shell script from /tree/bin:

```
generate_cram_csv.sh $crampath >> ${prefix}_cram.csv
```

[CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT](https://github.com/sanger-tol/treeval/blob/dev/modules/local/cram_filter_align_bwamem2_fixmate_sort.nf)
runs cram_filter and some conventional tools:

```
 // Please be aware one of the tools here required mem = 28 * reference size!!!
    """
    cram_filter -n ${from}-${to} ${cramfile} - | \\
        samtools fastq ${args1} | \\
        bwa-mem2 mem -p ${bwaprefix} -t${task.cpus} -5SPCp -H'${rglines}' - | \\
        samtools fixmate ${args3} - - | \\
        samtools sort ${args4} -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mem.bam -
```
[SAMTOOLS_MERGE](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/samtools/merge/main.nf) runs a samtools command:

```
    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${reference} \\
        ${prefix}.${file_type} \\
        $input_files
```

[pretext_map](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/pretextmap/main.nf) uses
*bioconda::pretextmap=0.1.9 bioconda::samtools=1.17* and logic depends on single or paired data:

```
 if [[ $input == *.pairs.gz ]]; then
        zcat $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    else
        samtools \\
            view \\
            $reference \\
            -h \\
            $input | PretextMap \\
            $args \\
            -o ${prefix}.pretext
    fi
```

[PRETEXT_INGEST_SNDRD](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/pretext_ingestion.nf) uses
[pretext_graph](https://github.com/sanger-tol/treeval/blob/dev/modules/local/pretext_graph.nf) with *conda "bioconda::pretextgraph=0.0.6 bioconda::ucsc-bigwigtobedgraph=448"* and runs:

```
//
// MODULE: PRETEXT GRAPH INGESTS THE OTHER TWO FILES DIRECTLY INTO THE PRETEXT
//          RUNNING AS IT'S OWN SUB IN ORDER TO NOT SLOW DOWN HIC_MAPPING ANY FURTHER
//

PRETEXT_GRAPH (
    pretext_file,
    ch_gap,
    coverage_file,
    cov_log_file,
    ch_telomere,
    repeat_cov_file
```

**PRETEXTMAP_HIGHRES** is a variant of pretextmap but called for high resolution, like the next optional **PRETEXT_INGEST_HIRES** step. These call the same command presumably with some changes to the parameters compared to the call above.

[SNAPSHOT_SRES](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/pretextsnapshot/main.nf) really uses the NF-core pretextsnapshot to make an image using:
```
 PretextSnapshot \\
        $args \\
        --map $pretext_map \\
        --prefix $prefix \\
        --folder .
```
[SAMTOOLS_MARKDUP](https://github.com/sanger-tol/treeval/tree/dev/modules/nf-core/samtools/markdup) calls samtools so probably use the existing samtools tool:
```
 samtools \\
        markdup \\
        $args \\
        ${reference} \\
        -@ $task.cpus \\
        -T $prefix \\
        $input \\
        ${prefix}.${extension}
```

[BAMTOBED_SORT](https://github.com/sanger-tol/treeval/blob/dev/modules/local/bamtobed_sort.nf) just uses samtools again to call:

```
samtools view -@${st_cores} -u -F0x400 ${bam} | bamToBed | sort -k4 --parallel=${task.cpus} -S ${buffer_mem}G > ${prefix}_merged_sorted.bed
```

[GET_PAIRED_CONTACT_BED](https://github.com/sanger-tol/treeval/blob/dev/modules/local/get_paired_contact_bed.nf) runs a /tree/bin shell script:

```
bed_to_contacts.sh $file > pre.bed
```

That script is some awk and sorting:

```
#!/bin/bash
paste -d '\t' - - < $1 | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else {print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '01'  | sort -k3,3d -k7,7d | awk 'NF==11'
```

If running the full workflow, there is a step using a /tree/bin java file that can be ignored for now
[JUICER_TOOLS_PRE](https://github.com/sanger-tol/treeval/blob/dev/modules/local/juicer_tools_pre.nf) calls:

```
java ${juicer_jvm_params} \\
        -jar ${projectDir}/bin/${juicer_tools_jar} pre \\
        ${pairs} \\
        ${prefix}.hic \\
        ${sizes}
```

[COOLER_CLOAD](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/cooler/cload/main.nf) uses *"bioconda::cooler=0.9.2"* and many cooler functions are available in the toolshed. It calls:

```
 cooler cload \\
        $args \\
        $nproc \\
        ${chromsizes}:${cool_bin} \\
        $pairs \\
        ${prefix}.cool
```
[COOLER_ZOOMIFY](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/cooler/zoomify/main.nf) uses the same dependency to run:

```
cooler zoomify \\
        $args \\
        -n $task.cpus \\
        -o ${prefix}.mcool \\
        $cool
```

