### [#6 generate_genome](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/generate_genome.nf)

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_generate_genome.png)

This executes three lines of shell script in [custom_get_chromsizes](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/custom/getchromsizes/main.nf) so probably needs a new tool:

```
ln -s ${fasta} ${prefix}.fa
    samtools faidx ${prefix}.fa -o ${prefix}.fa.fai
    cut -f 1,2 ${prefix}.fa.fai > ${prefix}.${suffix}
```

There’s a sort and the largest scaffold is found.

This can all be done with a new tool easily but is this available in any existing assembly workflows?

The scaffold size is related to Tabix - are there better ways of solving this?


```
workflow GENERATE_GENOME {
    take:
    reference_file  // Channel: path(file)

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    CUSTOM_GETCHROMSIZES (
        reference_file,
        "temp.genome"
    )
    ch_versions     = ch_versions.mix(  CUSTOM_GETCHROMSIZES.out.versions )

    //
    // MODULE: SORT CHROM SIZES BY CHOM SIZE NOT NAME
    //
    GNU_SORT (
        CUSTOM_GETCHROMSIZES.out.sizes
    )
    //
    // MODULE: Cut out the largest scaffold size and use as comparator against 512MB
    //          This is the cut off for TABIX using tbi indexes
    //
    GET_LARGEST_SCAFF (
        CUSTOM_GETCHROMSIZES.out.sizes
    )
```

### Potential tools to use:

* https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/fasta_compute_length/fasta_compute_length/1.0.3


