## [longread_coverage](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/longread_coverage.nf)

### Galaxy proposed solution

#### Partial workflow - missing log coverage and punchlists - TBA:
![Partial workflow](https://github.com/fubar2/treeval_gal/blob/main/longread_coverage/treevalgal_longread_coverage_wf0.png)

#### Sample JBrowse:

![Partial workflow](https://github.com/fubar2/treeval_gal/blob/main/longread_coverage/treevalgal_long_coverage_sample.png)

#### Zoomed in at a certian point gives a kind of pileup instead of a coverage histogram:

![Partial workflow](https://github.com/fubar2/treeval_gal/blob/main/longread_coverage/treevalgal_longread_coverage_zoom.png)



Note: bedtools does [coverage](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) so no need for the TreeVal perl script!

### Treeval NF DDL subworkflow deconstruction and explanation

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_longread_coverage.png)

Lots of DDL steps
Outputs are:
```
    treeval_upload/
        coverage.bw: Coverage of aligned reads across the reference genome in bigwig format.
        coverage_log.bw: A log corrected coverage file which aims to smooth out the above track.
    treeval_upload/punchlists/
        maxdepth.bigbed: Max read depth punchlist in bigBed format.
        zerodepth.bigbed: Zero read depth punchlist in bigBed format.
        halfcoverage.bigbed: Half read depth punchlist in bigBed format.
```

The DDL has function calls explained below.
Most of the rest of the DDL is not going to be needed other than to
figure out exactly how each function gets parameters supplied to the actual command lines.

```
workflow LONGREAD_COVERAGE {

    take:
    reference_tuple     // Channel: tuple [ val(meta), file( reference_file ) ]
    dot_genome          // Channel: tuple [ val(meta), [ file( datafile ) ]   ]
    reads_path          // Channel: tuple [ val(meta), val( str )             ]

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: CREATES INDEX OF REFERENCE FILE
    //
    MINIMAP2_INDEX(
        reference_tuple
    )
    ch_versions             = ch_versions.mix( MINIMAP2_INDEX.out.versions )

    //
    // LOGIC: PREPARE GET_READS_FROM_DIRECTORY INPUT
    //
    reference_tuple
        .combine( reads_path )
        .map { meta, ref, meta2, reads_path ->
            tuple(
                [   id          : meta.id,
                    single_end  : true  ],
                reads_path
            )
        }
        .set { get_reads_input }

    //
    // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
    //
    ch_grabbed_read_paths   = GrabFiles( get_reads_input )

    //
    // LOGIC: PACBIO READS FILES TO CHANNEL
    //
    ch_grabbed_read_paths
        .map { meta, files ->
            tuple( files )
        }
        .flatten()
        .set { ch_read_paths }

    //
    // LOGIC: COMBINE PACBIO READ PATHS WITH MINIMAP2_INDEX OUTPUT
    //
    MINIMAP2_INDEX.out.index
        .combine( ch_read_paths )
        .combine( reference_tuple )
        .map { meta, ref_mmi, read_path, ref_meta, reference ->
            tuple(
                [   id          : meta.id,
                    single_end  : true,
                    split_prefix: read_path.toString().split('/')[-1].split('.fasta.gz')[0]
                ],
                read_path,
                ref_mmi,
                true,
                false,
                false,
                file( reference ).size()
            )
        }
        .branch {
            large               : it[6] > 3000000000
            small               : it[6] < 3000000000
        }
        .set { mma_input }

    mma_input.large
        .multiMap { meta, read_path, ref_mmi, bam_output, cigar_paf, cigar_bam, file_size ->
            read_tuple          : tuple( meta, read_path)
            mmi_index           : ref_mmi
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { large }

    mma_input.small
        .multiMap { meta, read_path, ref_mmi, bam_output, cigar_paf, cigar_bam, file_size ->
            read_tuple          : tuple( meta, read_path)
            mmi_index           : ref_mmi
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { small }

    //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE <5GB PER SCAFFOLD
    //
    MINIMAP2_ALIGN (
        small.read_tuple,
        small.mmi_index,
        small.bool_bam_ouput,
        small.bool_cigar_paf,
        small.bool_cigar_bam
    )
    ch_versions             = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_align_bams           = MINIMAP2_ALIGN.out.bam

    //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE >5GB PER SCAFFOLD
    //
    MINIMAP2_ALIGN_SPLIT (
        large.read_tuple,
        large.mmi_index,
        large.bool_bam_ouput,
        large.bool_cigar_paf,
        large.bool_cigar_bam
    )
    ch_versions             = ch_versions.mix(MINIMAP2_ALIGN_SPLIT.out.versions)

    //
    // LOGIC: COLLECT OUTPUTTED BAM FILES FROM BOTH PROCESSES
    //
    ch_align_bams
        .mix( MINIMAP2_ALIGN_SPLIT.out.bam )
        .set { ch_bams }

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    ch_bams
        .map { meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [ id    : file[0].toString().split('/')[-1].split('_')[0] ], // Change sample ID
                file
            )
        }
        .set { collected_files_for_merge }

    //
    // MODULE: MERGES THE BAM FILES IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM
    SAMTOOLS_MERGE(
        collected_files_for_merge,
        reference_tuple,
        [[],[]]
    )
    ch_versions             = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // MODULE: SORT THE MERGED BAM BEFORE CONVERSION
    //
    SAMTOOLS_SORT (
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions             = ch_versions.mix( SAMTOOLS_MERGE.out.versions )

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    SAMTOOLS_SORT.out.bam
        .combine( reference_tuple )
        .multiMap { meta, bam, ref_meta, ref ->
                bam_input       :   tuple(
                                        [   id          : meta.id,
                                            sz          : bam.size(),
                                            single_end  : true  ],
                                        bam,
                                        []   // As we aren't using an index file here
                                    )
                ref_input       :   tuple(
                                        ref_meta,
                                        ref
                                    )
        }
        .set { view_input }
    //
    // MODULE: EXTRACT READS FOR PRIMARY ASSEMBLY
    //
    SAMTOOLS_VIEW(
        view_input.bam_input,
        view_input.ref_input,
        []
    )
    ch_versions             = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    //
    // MODULE: BAM TO PRIMARY BED
    //
    BEDTOOLS_BAMTOBED(
        SAMTOOLS_VIEW.out.bam
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    //
    // LOGIC: PREPARING Genome2Cov INPUT
    //
    BEDTOOLS_BAMTOBED.out.bed
        .combine( dot_genome )
        .multiMap { meta, file, my_genome_meta, my_genome ->
            input_tuple         :   tuple (
                                        [   id          :   meta.id,
                                            single_end  :   true    ],
                                        file,
                                        1
                                    )
            dot_genome          :   my_genome
            file_suffix         :   'bed'
        }
        .set { genomecov_input }

    //
    // MODULE: Genome2Cov
    //
    BEDTOOLS_GENOMECOV(
        genomecov_input.input_tuple,
        genomecov_input.dot_genome,
        genomecov_input.file_suffix
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    //
    // MODULE: SORT THE PRIMARY BED FILE
    //
    GNU_SORT(
        BEDTOOLS_GENOMECOV.out.genomecov
    )
    ch_versions             = ch_versions.mix(GNU_SORT.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    GETMINMAXPUNCHES(
        GNU_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix(GETMINMAXPUNCHES.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MAX(
        GETMINMAXPUNCHES.out.max
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_MERGE_MAX.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MIN(
        GETMINMAXPUNCHES.out.min
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_MERGE_MIN.out.versions)

    //
    // MODULE: GENERATE DEPTHGRAPH
    //
    GRAPHOVERALLCOVERAGE(
        GNU_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix(GRAPHOVERALLCOVERAGE.out.versions)
    ch_depthgraph           = GRAPHOVERALLCOVERAGE.out.part

    //
    // LOGIC: PREPARING FINDHALFCOVERAGE INPUT
    //
    GNU_SORT.out.sorted
        .combine( GRAPHOVERALLCOVERAGE.out.part )
        .combine( dot_genome )
        .multiMap { meta, file, meta_depthgraph, depthgraph, meta_my_genome, my_genome ->
            halfcov_bed     :       tuple( [ id : meta.id, single_end : true  ], file )
            genome_file     :       my_genome
            depthgraph_file :       depthgraph
        }
        .set { halfcov_input }

    //
    // MODULE: FIND REGIONS OF HALF COVERAGE
    //
    FINDHALFCOVERAGE(
        halfcov_input.halfcov_bed,
        halfcov_input.genome_file,
        halfcov_input.depthgraph_file
    )
    ch_versions             = ch_versions.mix(FINDHALFCOVERAGE.out.versions)

    //
    // LOGIC: PREPARING NORMAL COVERAGE INPUT
    //
    GNU_SORT.out.sorted
        .combine( dot_genome )
        .combine(reference_tuple)
        .multiMap { meta, file, meta_my_genome, my_genome, ref_meta, ref ->
            ch_coverage_bed :   tuple ([ id: ref_meta.id, single_end: true], file)
            genome_file     :   my_genome
        }
        .set { bed2bw_normal_input }

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG FOR NORMAL COVERAGE
    //
    BED2BW_NORMAL(
        bed2bw_normal_input.ch_coverage_bed,
        bed2bw_normal_input.genome_file
    )
    ch_versions             = ch_versions.mix(BED2BW_NORMAL.out.versions)

    //
    // MODULE: CONVERT COVERAGE TO LOG
    //
    LONGREADCOVERAGESCALELOG(
        GNU_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix(LONGREADCOVERAGESCALELOG.out.versions)

    //
    // LOGIC: PREPARING LOG COVERAGE INPUT
    //
    LONGREADCOVERAGESCALELOG.out.bed
        .combine( dot_genome )
        .combine(reference_tuple)
        .multiMap { meta, file, meta_my_genome, my_genome, ref_meta, ref ->
            ch_coverage_bed :   tuple ([ id: ref_meta.id, single_end: true], file)
            genome_file     :   my_genome
        }
        .set { bed2bw_log_input }

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG FOR LOG COVERAGE
    //
    BED2BW_LOG(
        bed2bw_log_input.ch_coverage_bed,
        bed2bw_log_input.genome_file
    )
    ch_versions             = ch_versions.mix(BED2BW_LOG.out.versions)

    //
    // LOGIC: GENERATE A SUMMARY TUPLE FOR OUTPUT
    //
    ch_grabbed_read_paths
            .collect()
            .map { meta, fasta ->
                tuple( [    id: 'pacbio',
                            sz: fasta instanceof ArrayList ? fasta.collect { it.size()} : fasta.size() ],
                            fasta
                )
            }
            .set { ch_reporting_pacbio }

```

[MINIMAP2_INDEX](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/minimap2/index/main.nf) uses *conda "bioconda::minimap2=2.24"* to execute

```
minimap2 \\
        -t $task.cpus \\
        -d ${fasta.baseName}.mmi \\
        $args \\
        $fasta
```
There is some DDL to split very large inputs before calling [MINIMAP2_ALIGN](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/minimap2/align/main.nf)
```
  //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE <5GB PER SCAFFOLD
    //
```

that uses *conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"*. It is called with "small" inputs to execute:

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
then the same MINIMAP2 function renamed as MINIMAP2_ALIGN_SPLIT is called with the large ones. No idea why this renamed function and logic to
adjustments ids and merge - there is a comment:
```
 //
    // MODULE: MERGES THE BAM FILES IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM
```

using [SAMTOOLS_MERGE](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/samtools/merge/main.nf) is called - also in [nuc_alignments](nuc_alignments) and [hic_mapping](hic_mapping)
uses *conda "bioconda::samtools=1.17"* to run

```
samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${reference} \\
        ${prefix}.${file_type} \\
        $input_files
```
[SAMTOOLS_SORT](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/samtools/sort/main.nf) follows using *conda "bioconda::samtools=1.17"* to execute:

```
 samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam
```

[SAMTOOLS_VIEW](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/samtools/view/main.nf) runs:

```
samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        ${readnames} \\
        $args \\
        -o ${prefix}.${file_type} \\
        $input \\
        $args2
```

[BEDTOOLS_BAMTOBED](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/bamtobed/main.nf) uses *conda "bioconda::bedtools=2.31.0"* to run:

```
 bedtools \\
        bamtobed \\
        $args \\
        -i $bam \\
        > ${prefix}.bed
```

[BEDTOOLS_GENOMECOV](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/genomecov/main.nf) uses *bioconda::bedtools=2.31.0"*.
That runs one of two separate command lines depending on whether the inputs are .bam files using odd DDL logic:

```
 if (intervals.name =~ /\.bam/) {
        """
        bedtools \\
            genomecov \\
            -ibam $intervals \\
            $args \\
            > ${prefix}.${extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        bedtools \\
            genomecov \\
            -i $intervals \\
            -g $sizes \\
            $args \\
            > ${prefix}.${extension}
```
[GNU_SORT](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/gnu/sort/main.nf) uses coreutils sort on the output merge called as:
```
 sort ${args} ${input} > ${output_file}
```

The sorted file is passed to [GETMINMAXPUNCHES](https://github.com/sanger-tol/treeval/blob/dev/modules/local/getminmaxpunches.nf) to run awk:

```
  $/
    cat "${bedfile}" \
    | awk '{ if ($4 == 0) {print $0 >> "zero.bed" } else if ($4 > 1000) {print $0 >> "max.bed"}}'
```
BEDTOOLS_MERGE_MIN and MAX are aliases for [BEDTOOLS_MERGE](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/bedtools/merge/main.nf) using *"bioconda::bedtools=2.31.0"* to execute:
```
bedtools \\
        merge \\
        -i $bed \\
        $args \\
        > ${prefix}.bed
```
It is called separately on the min/max output files from the GETMINMAXPUNCHES step above.

A depth graph is generated using [GRAPHOVERALLCOVERAGE](https://github.com/sanger-tol/treeval/blob/dev/modules/local/graphoverallcoverage.nf) using *"conda-forge::perl=5.26.2"* to run
a perl script from the tree/bin directory

```
graph_overall_coverage.pl $bed > ${prefix}.part
```

The script is:

```
#!/usr/bin/env perl

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

use warnings;

# my $file = shift;

my ($file) = @ARGV;

if (!@ARGV || ($ARGV[0] eq '--version')) {
    print "1.0\n";
    exit 0;
}

open (FILE, $file) || die "can't open file $file\n";

my %depthcount;
while (my $line = <FILE>) {
    chomp $line;
    my ($id, $start, $end, $depth) = split ("\t", $line);
    my $length = $end - $start;

    if ($depthcount{$depth}){
        $depthcount{$depth} += $length;
    }
    else {
        $depthcount{$depth} = $length;
    }
}

foreach my $depth (sort {$a<=>$b} keys %depthcount){
    print join("\t", $depth, $depthcount{$depth}) ."\n";
}
```
[FINDHALFCOVERAGE](https://github.com/sanger-tol/treeval/blob/dev/modules/local/findhalfcoverage.nf) runs a [python script](https://github.com/sanger-tol/treeval/blob/dev/bin/findHalfcoverage.py) from the bin directory
using the command:
```
  findHalfcoverage.py -c $bedfile -m $my_genome -d $depthgraph > ${prefix}.bed
```

GNU_SORT is used again, followed by BED2BW_NORMAL - that's an alias for [ucsc_bedgraphtobigwig](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/ucsc/bedgraphtobigwig/main.nf) that calls

```
    bedGraphToBigWig \\
            $bedgraph \\
            $sizes \\
            ${prefix}.bigWig
```

There is an existing bedgraph2BigWig wrapper: https://usegalaxy.eu/root?tool_id=wig_to_bigWig

Then there is a repeat of the above steps after log scaling with

[LONGREADCOVERAGESCALELOG](https://github.com/sanger-tol/treeval/blob/dev/modules/local/longreadcoveragescalelog.nf) that runs a python script from the /tree/bin directory using

BG: This script can be replaced by https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.0

```
longread_cov_log.py -i $bedfile > ${prefix}.bed
```
That python script is
```
#!/usr/bin/env python

import optparse
import math

# Script originally developed by Will Eagles (we3@sanger.ac.uk)


def process_line(line):
    line_values = line.rsplit(None, 1)

    try:
        cov_val = float(line_values[1])
    except:
        cov_val = 0

    if cov_val > 0:
        log_cov_val = math.log(cov_val)
    else:
        log_cov_val = 0

    return line_values[0] + "\t" + str(round(log_cov_val, 2))


def main():
    parser = optparse.OptionParser(version="%prog 1.0")
    parser.add_option(
        "-i",
        "--inputfile",
        dest="inputfile",
        default="default.input",
    )

    options, remainder = parser.parse_args()

    cov_bed = open(options.inputfile, "r")

    for line in cov_bed:
        print(process_line(line))


if __name__ == "__main__":
    main()
```

The output is processed with BED2BW_log, another alias for [ucsc_bedgraphtobigwig](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/ucsc/bedgraphtobigwig/main.nf) that calls

```
    bedGraphToBigWig \\
            $bedgraph \\
            $sizes \\
            ${prefix}.bigWig
```

