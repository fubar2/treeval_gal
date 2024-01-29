[PUNCHLIST](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/punchlist.nf)

Called by [nuc_alignments](nuc_alignments)

#### PAF is an optional output from minimap. Galaxy's minimap can produce a bam transparently to convert to a bed for viewing without pfaffing about with paf

This 33 line DDL subworkflow to make bed files from paf is largely redundant for Galaxy use.

```
#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { PAFTOOLS_SAM2PAF      } from '../../modules/nf-core/paftools/sam2paf/main'
include { PAF2BED               } from '../../modules/local/paf_to_bed'

workflow PUNCHLIST {
    take:
    reference_tuple // Channel: tuple [ val(meta), path(reference)]
    merged_bam      // Channel: tuple [ val(meta), path(bam_file)]

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: CONVERTS BAM INTO PAF FOR THE PUNCHLIST GENERATION
    //
    PAFTOOLS_SAM2PAF (
        merged_bam
    )
    ch_versions         = ch_versions.mix( PAFTOOLS_SAM2PAF.out.versions )

    //
    // MODULE: GENERATES PUNCHLIST FROM PAF FILE
    //
    PAF2BED (
        PAFTOOLS_SAM2PAF.out.paf
    )
    ch_versions         = ch_versions.mix( PAF2BED.out.versions )

    emit:
    punchlist           = PAF2BED.out.punchlist
    versions            = ch_versions.ifEmpty(null)
}
```

[PAFTOOLS_SAM2PAF](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/paftools/sam2paf/main.nf) uses samtools and
a minimap suite js ? 
```
samtools view -h ${bam} | paftools.js sam2paf - > ${prefix}.paf
```

[PAF2BED](https://github.com/sanger-tol/treeval/blob/dev/modules/local/paf_to_bed.nf) is a local module that calls the following command line :


```
paf_to_bed.sh ${file} ${prefix}.bed
```

That bash script is in /treeval/bin so a new tool is needed. Simply
```
#!/bin/bash

# paf_to_bed12.sh
# -------------------
# A shell script to convert a
# paf into bed format for use
# in JBrowse
# -------------------
# Author = yy5

version='1.0.0'

if [ $1 == '-v'];
then
    echo "$version"
else
    cat $1 | awk 'BEGIN{FS="\t";}{a[$1]++;if(a[$1]==2)print v[$1] ORS $0;if(a[$1]>2)print;v[$1]=$0;}' | awk '$(NF+1) = ($10/$11)*100' | awk '$(NF+1) = ($10/$2)*100' | awk -vOFS='\t' '{print $6,$8,$9,$1,$2,$10,$(NF-1),$NF}' > $2
fi
```
