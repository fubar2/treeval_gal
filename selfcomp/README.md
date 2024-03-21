# treeval_gal
Sanger treeval nf workflow translation into Galaxy work in progress

Using mummer4 allows larger contigs but the smaller contigs are apparently a good idea.
Without them the mapped bams might be 5GB so too big for JB2 but a bed or bigwig is useable.
Yumi Sims suggests that the homology mapping is improved by splitting the fasta into smaller segments.
Anyone know why?

Not sure if we need to block stitching step below - could be done with bedtools if needed...

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_self_comp.png)


```
The selfcomp subworkflow is a comparative genomics analysis algorithm originally performed by the Ensembl projects database, and reverse engineered in Python3 by @yumisims. It involves comparing the genes and genomic sequences within a single species. The goal of the analysis is to identify haplotypic duplications in a particular genome assembly.
Output files

    treeval_upload/
        *_selfcomp.bigBed: BigBed file containing selfcomp track data.

```

```
#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { MUMMER                         } from '../../modules/nf-core/mummer/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UCSC_BEDTOBIGBED               } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                  } from '../../modules/nf-core/bedtools/sort/main'
include { SELFCOMP_SPLITFASTA            } from '../../modules/local/selfcomp_splitfasta'
include { SELFCOMP_MUMMER2BED            } from '../../modules/local/selfcomp_mummer2bed'
include { SELFCOMP_MAPIDS                } from '../../modules/local/selfcomp_mapids'
include { CHUNKFASTA                     } from '../../modules/local/chunkfasta'
include { CONCATMUMMER                   } from '../../modules/local/concatmummer'
include { SELFCOMP_ALIGNMENTBLOCKS       } from '../../modules/local/selfcomp_alignmentblocks'
include { CONCATBLOCKS                   } from '../../modules/local/concatblocks'
include { BEDTOOLS_MERGE                 } from '../../modules/nf-core/bedtools/merge/main'

workflow SELFCOMP {
    take:
    reference_tuple      // Channel: tuple [ val(meta), path(reference_file) ]
    dot_genome           // Channel: tuple [ val(meta), [ path(datafile) ] ]
    mummer_chunk         // Channel: val( int )
    motif_len            // Channel: val( int )
    selfcomp_as          // Channel: val( dot_as location )

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: SPLITS INPUT FASTA INTO 500KB CHUNKS
    //         EMITS CHUNKED FASTA
    //
    SELFCOMP_SPLITFASTA(
        reference_tuple
    )
    ch_versions             = ch_versions.mix( SELFCOMP_SPLITFASTA.out.versions )

    //
    // MODULE: SPLIT INPUT FASTA INTO 1GB CHUNKS
    //         EMITS CHUNKED FASTA
    //
```

The input assembly is split with a perl script in [selfcomp split fasta](https://github.com/sanger-tol/treeval/blob/dev/bin/split_genomes_for_ensembl.pl)

```
#!/usr/bin/env perl

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

#
# Take a genome fasta file, and creates a fasta file and agp
# that is suitable for loading into Ensembl
#
# If the genome is already suitable for loading (i.e. all of the
# sequences are small enough for direct loading), then the script
# says so and does nothing
#

use strict;

use Bio::SeqIO;
use Getopt::Long;

my $MAX_CONTIG_LEN = 500000;

my ($genome_fa,
    $contig_fa_file,
    $agp_file) = @ARGV;

if (!@ARGV || ($ARGV[0] eq '--version')) {
    print "1.0\n";
    exit 0;
}

my $seqio = Bio::SeqIO->new(-format  => 'fasta',
                            -file => ($genome_fa =~ /\.gz$/) ? "gunzip -c $genome_fa |" : $genome_fa,
);

my (@toplevels, $need_agp, $count);

while (my $seq = $seqio->next_seq) {
    if ($seq->length > $MAX_CONTIG_LEN) {
        $need_agp = 1;
    }
    push @toplevels, $seq;
}

if (not $need_agp) {
    print "All sequences are small enough for direct loading. Did not create AGP\n";
} else {
    my $outseqio = Bio::SeqIO->new(-format => 'fasta',
                                    -file => ">$contig_fa_file");
    open(my $agp_fh, ">$agp_file") or die "Could not open $agp_file for writing\n";

    foreach my $seq (@toplevels) {
        if ($seq->length < $MAX_CONTIG_LEN) {
            $outseqio->write_seq($seq);
            printf($agp_fh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n", $seq->id, 1, $seq->length, ++$count, $seq->id, 1, $seq->length);
        } else {
            print STDERR "GOT ONE\n";
            my $seg_count = 1;
            for(my $i = 0; $i < $seq->length; $i += $MAX_CONTIG_LEN) {
                my $start = $i + 1;
                my $end = $start + $MAX_CONTIG_LEN - 1;
                $end = $seq->length if $end > $seq->length;

                my $new_id = sprintf("%s.%d", $seq->id, $seg_count++);

                my $seq_seg = Bio::PrimarySeq->new( -id => $new_id,
                                                    -seq => substr($seq->seq, $start - 1, $end - $start + 1));
                $outseqio->write_seq($seq_seg);

                printf($agp_fh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n", $seq->id, $start, $end, ++$count, $seq_seg->id, 1, $seq_seg->length);
            }
        }
    }
}

exit(0);
```

```
    CHUNKFASTA(
        SELFCOMP_SPLITFASTA.out.fa,
        mummer_chunk
    )
    ch_versions             = ch_versions.mix( CHUNKFASTA.out.versions )

```


Fasta is [chunked](https://github.com/sanger-tol/treeval/blob/dev/modules/local/chunkfasta.nf) with conda pyfasta=0.5.2-1

```
    script:
    """
    pyfasta split -n $number_of_chunks $fasta

```

Then it's mapped with mummer.


```
    //
    // LOGIC: CONVERTS ABOVE OUTPUTS INTO A SINGLE TUPLE
    //
    ch_query_tup = CHUNKFASTA.out.fas
        .map{ meta, query ->
            [query]
        }
        .flatten()

    ch_ref = SELFCOMP_SPLITFASTA.out.fa
        .map{ meta, ref ->
            ref
        }

    ch_mummer_input = ch_query_tup
        .combine(ch_ref)
        .map{ query, ref ->
                tuple([   id: query.toString().split('/')[-1] ],
                        ref,
                        query
                )
        }

    //
    // MODULE: ALIGNS 1GB CHUNKS TO 500KB CHUNKS
    //         EMITS MUMMER ALIGNMENT FILE
    //
    MUMMER(
        ch_mummer_input
    )
    ch_versions             = ch_versions.mix( MUMMER.out.versions )

    //
    // LOGIC: GROUPS OUTPUT INTO SINGLE TUPLE BASED ON REFERENCE META
    //
    MUMMER.out.coords
        .combine( reference_tuple )
        .map { coords_meta, coords, ref_meta, ref ->
                tuple(  ref_meta,
                        coords
                )
        }
        .groupTuple( by:[0] )
        .set{ ch_mummer_files }


    //
    // MODULE: MERGES MUMMER ALIGNMENT FILES
    //
    CONCATMUMMER(
        ch_mummer_files
    )
    ch_versions             = ch_versions.mix( CONCATMUMMER.out.versions )
```

Mummer `delta` outputs are used because they have Mummer 3 - nucmer in Galaxy now has Mummer 4 with bam/cram outputs.
They're text [so cat is used](https://github.com/sanger-tol/treeval/blob/dev/modules/local/concatmummer.nf) to join them.

```
 cat $coords > ${prefix}.mummer
```
That `delta` is converted to `bed` and all the split contigs and positions are fixed up.


```
    //
    // MODULE: CONVERT THE MUMMER ALIGNMENTS INTO BED FORMAT
    //
    SELFCOMP_MUMMER2BED(
        CONCATMUMMER.out.mummer,
        motif_len
    )
    ch_versions             = ch_versions.mix( SELFCOMP_MUMMER2BED.out.versions )

    //
    // MODULE: GENERATE A LIST OF IDs AND GENOMIC POSITIONS OF SELFCOMPLEMENTARY REGIONS
    //         EMITS BED FILE
    //
    SELFCOMP_MAPIDS(
        SELFCOMP_MUMMER2BED.out.bedfile,
        SELFCOMP_SPLITFASTA.out.agp
    )
    ch_versions             = ch_versions.mix( SELFCOMP_MAPIDS.out.versions )

    //
    // MODULE: SORTS ABOVE OUTPUT BED FILE AND RETAINS BED SUFFIX
    //
    BEDTOOLS_SORT(
        SELFCOMP_MAPIDS.out.bedfile,
        []
    )
    ch_versions             = ch_versions.mix( BEDTOOLS_SORT.out.versions )

    //
    // MODULE: BUILD ALIGNMENT BLOCKS
    //
    SELFCOMP_ALIGNMENTBLOCKS(
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix( SELFCOMP_ALIGNMENTBLOCKS.out.versions )

    //
    // MODULE: SORT BLOCKS FILES AND FILTER BY MOTIF LENGTH
    //
    CONCATBLOCKS(
        SELFCOMP_ALIGNMENTBLOCKS.out.blockfile
    )
    ch_versions             = ch_versions.mix( CONCATBLOCKS.out.versions )
```

The [last stage](https://github.com/sanger-tol/treeval/blob/dev/bin/build_alignment_block.py) is a python script to stitch bed blocks using pandas.

```
#!/usr/bin/env python

import pandas as pd
import os
import sys
import string
import random
import csv
import optparse
import pybedtools
from pybedtools import BedTool


def sort_blocks(df):
    return df.sort_values(["rstart", "qstart"], ascending=[True, True])


def get_block(df, index):
    block = pd.DataFrame([])
    if index < (len(small_cluster_sort.index) - 2):
        if df.iloc[index].loc["qstart"] > df.iloc[index + 1].loc["qstart"]:
            block = df[0 : index + 1]
            leftover = df[index + 1 : len(df.index) - 1]
            qmin = leftover[["qstart"]].min()
            print(qmin)
            index_list = list(range(0, index + 1))
            df.drop(df.index[index_list], inplace=True)

    return block, df


def arrange_fields(df):
    df["fragid"] = df["qchr"].str.cat(df["qstart"].astype(str), sep=":").str.cat(df["qend"].astype(str), sep=":")

    return df[["refchr", "rstart", "rend", "fragid", "qstrand"]]


def build_block(mylist):
    qmin = 0
    qlist = []
    nlist = []

    for idx, x in enumerate(mylist):
        if idx < len(mylist) - 1:
            qcurrent = int(x[6])
            rcurrent = int(x[1])
            qnext = mylist[idx + 1][6]
            leftover = mylist[idx : len(mylist)]

            # leftd = int(max((x[6]-qcurrent for x in leftover), default=0))

            leftd = list(x[6] - qmin for x in leftover)

            positives = [x for x in leftd if x > 0]

            min_value = min((positives), default=0)

            indmin = leftd.index(min_value)

            rm = leftover[indmin][1]

            if qcurrent > qmin and qcurrent < qnext and rm == rcurrent:
                qmin = qcurrent
                qlist.append(idx)

            if qcurrent > qmin and qcurrent < qnext and rm > rcurrent:
                nlist.append(idx)

            if qcurrent > qmin and qcurrent > qnext:
                nlist.append(idx)

            if qcurrent < qmin and qcurrent > qnext:
                nlist.append(idx)

        if idx == len(mylist) - 1:
            if mylist[idx][6] > qmin:
                qlist.append(idx)
            else:
                nlist.append(idx)

    alignment_chain = [mylist[i] for i in qlist]
    new_list = [mylist[i] for i in nlist]
    return alignment_chain, new_list


#########main##########

parser = optparse.OptionParser(version="%prog 1.0")
parser.add_option(
    "-i",
    "--input",
    dest="input_bedfile",
    default="default.input",
)
parser.add_option(
    "-o",
    "--output",
    dest="out_bedfile",
    default="default.output",
)
options, remainder = parser.parse_args()

inbed = options.input_bedfile
outbed = options.out_bedfile

fo = open(outbed, "a")

sc = pd.read_csv(inbed, sep="\t", comment="#", header=None)

sc.columns = ["refchr", "rstart", "rend", "qchr", "maplen", "qstrand", "qstart", "qend", "qlen"]

ans = [y for x, y in sc.groupby("refchr")]

for mycluster in ans:
    for small_cluster in [y for x, y in mycluster.groupby("qchr")]:
        small_cluster_sort = sort_blocks(small_cluster)

        newdf = small_cluster_sort.reset_index(drop=True)

        newlist = newdf.values.tolist()

        while newlist:
            blocks, newlist = build_block(newlist)

            # fileprefix = "".join(random.choices(string.ascii_lowercase + string.digits, k=12))
            # filename = fileprefix + ".block"
            newblocks = [
                [x if i != 3 else y[3] + ":" + str(y[6]) + ":" + str(y[7]) for i, x in enumerate(y)] for y in blocks
            ]

            a = pybedtools.BedTool(newblocks)
            merged = a.merge(d=100000, c="4,7,8", o="collapse,min,max", delim="|")
            fo.write(str(merged))

fo.close()
```

```
    //
    // MODULE: CONVERTS ABOVE OUTPUT INTO BIGBED FORMAT
    //
    UCSC_BEDTOBIGBED(
        CONCATBLOCKS.out.chainfile,
        dot_genome.map{it[1]}, // Pulls file from tuple ( meta and file )
        selfcomp_as
    )
    ch_versions             = ch_versions.mix( UCSC_BEDTOBIGBED.out.versions )

    emit:
    ch_bigbed               = UCSC_BEDTOBIGBED.out.bigbed
    versions                = ch_versions.ifEmpty(null)
}
```
