
### [#4 gap_finder](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/gap_finder.nf)

This writes a new file in a DDL specific way so need some help confirming the logic.
![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_gap_finder.png)
```
workflow GAP_FINDER {
    take:
    reference_tuple     // Channel [ val(meta), path(fasta) ]
    max_scaff_size      // val(size of largest scaffold in bp)
    main:
    ch_versions     = Channel.empty()
    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        reference_tuple
    )
    ch_versions     = ch_versions.mix( SEQTK_CUTN.out.versions )
    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAP_LENGTH (
        SEQTK_CUTN.out.bed
    )
    ch_versions     = ch_versions.mix( GAP_LENGTH.out.versions )
    //
    // LOGIC: Adding the largest scaffold size to the meta data so it can be used in the modules.config
    //
    SEQTK_CUTN.out.bed
        .combine(max_scaff_size)
        .map {meta, row, scaff ->
            tuple([ id          : meta.id,
                    max_scaff   : scaff >= 500000000 ? 'csi': ''
                ],
                file(row)
            )}
        .set { modified_bed_ch }
    //
    // MODULE: BGZIP AND TABIX THE GAP FILE
    //
    TABIX_BGZIPTABIX (
        modified_bed_ch
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    gap_file        = GAP_LENGTH.out.bedgraph
    gap_tabix       = TABIX_BGZIPTABIX.out.gz_csi
    versions        = ch_versions.ifEmpty(null)
}
```


[seqtk is available](https://toolshed.g2.bx.psu.edu/view/iuc/package_seqtk_1_0_r75/8b7d6f6cb89b) in the Toolshed. [SEQTK_CUTN](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/seqtk/cutn/main.nf) is a nf-core module and it just does this:


```
    seqtk \\
            cutN \\
            $args \\
            -g $fasta \\
            > ${prefix}.bed
```


[GAP_LENGTH](https://github.com/sanger-tol/treeval/blob/dev/modules/local/gap_length.nf) is a local nf module that runs awk. If _a = $3-$2_ the math seems to return sqrt(a*a) or abs(a) - awk does not have an abs function, to enforce positive values, on the output from seqtk.:


```
    $/
    cat "${file}" \
    | awk '{print $0"\t"sqrt(($3-$2)*($3-$2))}' > ${prefix}_gap.bedgraph
```


Should not be hard to implement.

There is some intriguing DDL that seems to write a file with metadata added to each row of an intermediate bed file. Need some content expertise to understand, duplicate and test this.


```
    SEQTK_CUTN.out.bed
            .combine(max_scaff_size)
            .map {meta, row, scaff ->
                tuple([ id          : meta.id,
                        max_scaff   : scaff >= 500000000 ? 'csi': ''
                    ],
                    file(row)
                )}
            .set { modified_bed_ch }
```
[TABIX_BGZIPTABIX](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/gap_finder.nf) uses *bioconda::tabix=1.11* preceded by 
some interesting DDL that seems to flag the need for csi indexes if largest scaffold is too big to use tabix:

```
 SEQTK_CUTN.out.bed
        .combine(max_scaff_size)
        .map {meta, row, scaff ->
            tuple([ id          : meta.id,
                    max_scaff   : scaff >= 500000000 ? 'csi': ''
                ],
                file(row)
            )}
        .set { modified_bed_ch }
    //
```
The actual call is just:

```
 bgzip  --threads ${task.cpus} -c $args $input > ${prefix}.${input.getExtension()}.gz
 tabix $args2 ${prefix}.${input.getExtension()}.gz
```

The 500MB max_scaff limit is a tabix limit apparently….
```
The tabix (.tbi) and BAI index formats can handle individual chromosomes up to 512 Mbp (2^29 bases) in length. If your input file might contain data lines with begin or end positions greater than that, you will need to use a CSI index.
```

Tabix bgzip [is available from the IUC](https://toolshed.g2.bx.psu.edu/repository/browse_repositories?f-free-text-search=tabix&sort=name&operation=view_or_manage_repository&id=2c71da8851968c89) in the Toolshed to run on that output file.
