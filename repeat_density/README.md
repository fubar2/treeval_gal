### [#15 repeat_density](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/repeat_density.nf)



![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_repeat_density.png)


```
workflow REPEAT_DENSITY {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    dot_genome
    main:
    ch_versions         = Channel.empty()
    //
    // MODULE: MARK UP THE REPEAT REGIONS OF THE REFERENCE GENOME
    //
    WINDOWMASKER_MKCOUNTS (
        reference_tuple
    )
    ch_versions         = ch_versions.mix( WINDOWMASKER_MKCOUNTS.out.versions )
    //
    // MODULE: CALCULATE THE STATISTICS OF THE MARKED UP REGIONS
    //
    WINDOWMASKER_USTAT(
        WINDOWMASKER_MKCOUNTS.out.counts,
        reference_tuple
    )
    ch_versions         = ch_versions.mix( WINDOWMASKER_USTAT.out.versions )
    //
    // MODULE: USE USTAT OUTPUT TO EXTRACT REPEATS FROM FASTA
    //
    EXTRACT_REPEAT(
        WINDOWMASKER_USTAT.out.intervals
    )
    ch_versions         = ch_versions.mix( EXTRACT_REPEAT.out.versions )
    //
    // MODULE: CREATE WINDOWS FROM .GENOME FILE
    //
    BEDTOOLS_MAKEWINDOWS(
        dot_genome
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_MAKEWINDOWS.out.versions )
    //
    // LOGIC: COMBINE TWO CHANNELS AND OUTPUT tuple(meta, windows_file, repeat_file)
    //
    BEDTOOLS_MAKEWINDOWS.out.bed
        .combine( EXTRACT_REPEAT.out.bed )
        .map{ meta, windows_file, repeat_meta, repeat_file ->
                    tuple (
                        meta,
                        windows_file,
                        repeat_file
                    )
        }
        .set { intervals }
    //
    // MODULE: GENERATES THE REPEAT FILE FROM THE WINDOW FILE AND GENOME FILE
    //
    BEDTOOLS_INTERSECT(
        intervals,
        dot_genome
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_INTERSECT.out.versions )
    //
    // MODULE: FIXES IDS FOR REPEATS
    //
    RENAME_IDS(
        BEDTOOLS_INTERSECT.out.intersect
    )
    ch_versions         = ch_versions.mix( RENAME_IDS.out.versions )
    //
    // MODULE: SORTS THE ABOVE BED FILES
    //
    GNU_SORT_A (
        RENAME_IDS.out.bed              // Intersect file
    )
    ch_versions         = ch_versions.mix( GNU_SORT_A.out.versions )
    GNU_SORT_B (
        dot_genome                      // Genome file - Will not run unless genome file is sorted to
    )
    ch_versions         = ch_versions.mix( GNU_SORT_B.out.versions )
    GNU_SORT_C (
        BEDTOOLS_MAKEWINDOWS.out.bed    // Windows file
    )
    ch_versions         = ch_versions.mix( GNU_SORT_C.out.versions )
    //
    // MODULE: ADDS 4TH COLUMN TO BED FILE USED IN THE REPEAT DENSITY GRAPH
    //
    REFORMAT_INTERSECT (
        GNU_SORT_A.out.sorted
    )
    ch_versions         = ch_versions.mix( GNU_SORT_C.out.versions )
    //
    // LOGIC: COMBINES THE REFORMATTED INTERSECT FILE AND WINDOWS FILE CHANNELS AND SORTS INTO
    //      tuple(intersect_meta, windows file, intersect file)
    //
    REFORMAT_INTERSECT.out.bed
        .combine( GNU_SORT_C.out.sorted )
        .map{ intersect_meta, bed, sorted_meta, windows_file ->
                    tuple (
                        intersect_meta,
                        windows_file,
                        bed
                    )
        }
        .set { for_mapping }
    //
    // MODULE: MAPS THE REPEATS AGAINST THE REFERENCE GENOME
    //
    BEDTOOLS_MAP(
        for_mapping,
        GNU_SORT_B.out.sorted
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_MAP.out.versions )
    //
    // MODULE: REPLACES . WITH 0 IN MAPPED FILE
    //
    REPLACE_DOTS (
        BEDTOOLS_MAP.out.map
    )
    ch_versions         = ch_versions.mix( REPLACE_DOTS.out.versions )
    //
    // MODULE: CONVERTS GENOME FILE AND BED INTO A BIGWIG FILE
    //
    UCSC_BEDGRAPHTOBIGWIG(
        REPLACE_DOTS.out.bed,
        GNU_SORT_B.out.sorted.map { it[1] } // Pulls file from tuple of meta and file
    )
    ch_versions         = ch_versions.mix( UCSC_BEDGRAPHTOBIGWIG.out.versions )
    emit:
    repeat_density      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    versions            = ch_versions.ifEmpty(null)
}
```


First step, [WINDOWMASKER_MKCOUNTS](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/windowmasker/mk_counts/main.nf)
uses an NF core module to run [a tool already available ](https://toolshed.g2.bx.psu.edu/view/yating-l/windowmasker_2_5_0/f80c9e6700ba) in the toolshed
with this command line:


```
     windowmasker -mk_counts \\
            $args \\
            -mem ${memory} \\
            -in ${ref} \\
            -out ${prefix}.txt
```


That repository also includes the command used in the next **[ustat](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/windowmasker/mk_counts/main.nf)** step.


```
     windowmasker -ustat \\
            ${counts} \\
            $args \\
            -in ${ref} \\
            -out ${output}
```

https://github.com/goeckslab/WindowMasker may be useful here.

[extract_repeats](https://github.com/sanger-tol/treeval/blob/dev/modules/local/extract_repeat.nf) runs a perl script in the /tree/bin directory


```
extract_repeat.pl $file > ${prefix}_repeats.bed
```


so a new simple tool is needed.

Then there are calls to [makewindows](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/makewindows/main.nf)

```
        bedtools \\
            makewindows \\
            ${arg_input} \\
            ${args} \\
            > ${prefix}.bed
```


and [intersect](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/intersect/main.nf)


```
        bedtools \\
            intersect \\
            -a $intervals1 \\
            -b $intervals2 \\
            $args \\
            $sizes \\
            > ${prefix}.${extension}
```


 that can probably be done with existing toolshed bedtools.

[rename_ids](https://github.com/sanger-tol/treeval/blob/dev/modules/local/rename_ids.nf) is a simple sed operation to replace a
period with 0 - replace_dots below does the same thing.


```
        $/
        cat "${file}" \
        | sed 's/\./0/g' > ${prefix}_renamed.bed
```


so could use the sed tool or a new specialised tool.

[reformat_intersect](https://github.com/sanger-tol/treeval/blob/dev/modules/local/reformat_intersect.nf) is a local NF module that does


```
     $/
        cat "${file}" \
        | awk '{print $0"\t"sqrt(($3-$2)*($3-$2))}'\
        | sed 's/\./0/g' > ${prefix}.bed
```


so a new tool perhaps. Odd that it does replace_dots sed - is it redundant below?

[Bedtools map](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/map/main.nf) does:


```
    bedtools \\
            map \\
            -a $intervals1 \\
            -b $intervals2 \\
            $args \\
            $sizes \\
            > ${prefix}.${extension}
```


so can use the existing bedtools tools.

[replace_dots ](https://github.com/sanger-tol/treeval/blob/dev/modules/local/replace_dots.nf) (might be redundant give reform_intersect above but) calls the following sed command:


```
     $/
        cat "${file}" | sed 's/\./0/g' > "${prefix}_nodot.bed"
```


so a new tool will be needed or the IUC sed tool perhaps.

[ucsc_bedgraphtobigwig](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/ucsc/bedgraphtobigwig/main.nf) calls


```
    bedGraphToBigWig \\
            $bedgraph \\
            $sizes \\
            ${prefix}.bigWig
```


We do have a bedgraph2BigWig wrapper: https://usegalaxy.eu/root?tool_id=wig_to_bigWig
~There is a [tool in the toolshed for ribogalaxy](https://toolshed.g2.bx.psu.edu/view/jackcurragh/ribogalaxy_bedgraphtobigwig/66f3a55ff2ec) that may be usable
or else a new tool will be needed.~

