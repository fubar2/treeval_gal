### [#2 ancestral_gene](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/ancestral_gene.nf)

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_ancestral_gene.png)

Calls 4 functions (modules or subworkflows) in **bold** below:


```
        //
    // MODULE: EXTRACTS ANCESTRALLY LINKED BUSCO GENES FROM FULL TABLE
    //
    EXTRACT_ANCESTRAL(
        ch_grab,
        ancestral_table
    )
    ch_versions             = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)
    //
    // LOGIC: STRIP OUT METADATA
    //
    ch_grab
        .map { meta, fulltable
                -> fulltable
            }
        .set { assignanc_input }
    //
    // MODULE: ASSIGN EXTRACTED GENES TO ANCESTRAL GROUPS
    //
    ASSIGN_ANCESTRAL(
        EXTRACT_ANCESTRAL.out.comp_location,
        assignanc_input
    )
    ch_versions             = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)
    //
    // MODULES: SORT THE BED FILE
    //
    BEDTOOLS_SORT(
        ASSIGN_ANCESTRAL.out.assigned_bed,
        []
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_SORT.out.versions)
    //
    // MODULES: CONVERT BED TO INDEXED BIGBED
    //
    UCSC_BEDTOBIGBED(
        BEDTOOLS_SORT.out.sorted,
        dot_genome.map{ it[1] },    // Pull file from tuple(meta, file)
        buscogene_as
    )
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)
    emit:
    ch_ancestral_bigbed     = UCSC_BEDTOBIGBED.out.bigbed
    versions                = ch_versions.ifEmpty(null)
```


The first step calls a local NF DDL module **[EXTRACT_ANCESTRAL](https://github.com/sanger-tol/treeval/blob/dev/modules/local/extract_ancestral.nf)** - mostly setup code in preparation for running the following command:

**_     buscopainter.py -r $ancestraltable -q $fulltable _**

[buscopainter ](https://github.com/charlottewright/buscopainter)is not part of [Busco,](https://toolshed.g2.bx.psu.edu/repository?repository_id=99f4c45aebf38997&changeset_revision=2a5b8b9936bf) so a new tool is needed with a content expert to advise and test. Need a content expert to help figure out how these parameters should be setup in Galaxy tool XML.

The second step, **[ASSIGN_ANCESTRAL](https://github.com/sanger-tol/treeval/blob/dev/modules/local/assign_ancestral.nf)**, is a local NF module. It runs a one line bash script:

**_assign_anc.py -l $comp_location -f $fulltable -c ${prefix}_assigned.bed_**

[That python code](https://github.com/sanger-tol/treeval/blob/f8e4b3bbbd75be6fa7ea6788337664d2533cdbdb/bin/assign_anc.py) is also found in the treeval/bin directory. Again, a new tool is needed and a content expert to advise and to test. The _$variables_ are DDL but work for Galaxy tools too._ ${prefix}_ is a subtask name idiom in DDL. Need a content expert to help figure out how these parameters should be setup in Galaxy tool XML.

The third step is a bedtools sort on that output, already available as an IUC bedtools tool.

The last step is [UCSC-bedtobigbed](https://open.bioqueue.org/home/knowledge/showKnowledge/sig/ucsc-bedtobigbed) using:

**_ `bedToBigBed \\`_**


```
        $bed \\
        $sizes \\
        $as_option \\
        $args \\
        ${prefix}.bigBed
```


Not found in the toolshed, so again, another trivial new tool needed and expert advice on how to get those parameters from users. It can probably be decoded by looking at all the DDL carefully but frankly, life is too short.

A quick check of the Toolshed shows that bedtools sort ([sort_bed.xml](https://toolshed.g2.bx.psu.edu/repository/browse_repository?id=8d84903cc667dbe7#)) is available, but 3 other requirements do not seem to be findable, so 3 new Galaxy tools need to be built and acceptance-tested by a content expert:



1. [buscopainter ](https://github.com/charlottewright/buscopainter)
2. [assign_ancestry](https://github.com/sanger-tol/treeval/blob/f8e4b3bbbd75be6fa7ea6788337664d2533cdbdb/bin/assign_anc.py)
3. [ucsc bedtobigbed](http://UCSC-bedtobigbed)

They are mostly one line scripts to be turned into simple new tools. Finding all the necessary test data will be more work.

