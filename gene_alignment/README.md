### [#5 gene_alignment](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/gene_alignment.nf)

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_gene_alignment.png)

This one is complicated and content expertise is needed to describe the logic properly. Looks like a few new tools too. I am just guessing at the DDL - the flow diagram helps.


```
    // This subworkflow takes an input fasta sequence and csv style list of organisms to return
    // bigbed files containing alignment data between the input fasta and csv style organism names.
    // Input - Assembled genomic fasta file
    // Output - A BigBed file per datatype per organism entered via csv style in the yaml.
    //
    // SUBWORKFLOW IMPORT BLOCK
    //
    include { PEP_ALIGNMENTS                    } from './pep_alignments'
    include { NUC_ALIGNMENTS as GEN_ALIGNMENTS  } from './nuc_alignments'
    include { NUC_ALIGNMENTS as RNA_ALIGNMENTS  } from './nuc_alignments'
    include { NUC_ALIGNMENTS as CDS_ALIGNMENTS  } from './nuc_alignments'
    workflow GENE_ALIGNMENT {
    take:
    dot_genome          // Channel [ val(meta), path(file) ]
    reference_tuple     // Channel [ val(meta), path(file) ]
    reference_index     // Channel [ val(meta), path(file) ]
    max_scaff_size      // Channel val(size of largest scaffold in bp)
    assembly_classT     // Channel val(clade_id)
    alignment_datadir   // Channel val(geneset_dir)
    alignment_genesets  // Channel val(geneset_id)
    alignment_common    // Channel val(common_name) // Not yet in use
    intron_size         // Channel val(50k)
    as_files            // Channel [ val(meta), path(file) ]
    main:
    ch_versions         = Channel.empty()
    //
    // LOGIC: TAKES A SINGLE LIKE CSV STRING AND CONVERTS TO LIST OF VALUES
    //          LIST IS MERGED WITH DATA_DIRECTORY AND ORGANISM_CLASS
    //
    ch_data             = alignment_genesets
                            .splitCsv()
                            .flatten()
    //
    // LOGIC:   COMBINE CH_DATA WITH ALIGNMENT_DIR AND ASSEMBLY_CLASS
    //          CONVERTS THESE VALUES INTO A PATH AND DOWNLOADS IT, THEN TURNS IT TO A TUPLE OF
    //          [ [ META.ID, META.TYPE, META.ORG ], GENE_ALIGNMENT_FILE ]
    //          DATA IS THEN BRANCHED BASED ON META.TYPE TO THE APPROPRIATE
    //          SUBWORKFLOW
    //
    ch_data
        .combine( alignment_datadir )
        .combine( assembly_classT )
        .map {
            ch_org, data_dir, classT ->
                file("${data_dir}${classT}/csv_data/${ch_org}-data.csv")
        }
        .splitCsv( header: true, sep:',')
        .map( row ->
        tuple([ org:    row.org,
                type:   row.type,
                id:     row.data_file.split('/')[-1].split('.MOD.')[0]
            ],
            file(row.data_file)
        ))
        .branch {
            pep: it[0].type  == 'pep'
            gen: it[0].type  == 'cdna'
            rna: it[0].type  == 'rna'
            cds: it[0].type  == 'cds'
        }
        .set {ch_alignment_data}
    pep_files = ch_alignment_data.pep.collect()
    gen_files = ch_alignment_data.gen.collect()
    rna_files = ch_alignment_data.rna.collect()
    cds_files = ch_alignment_data.cds.collect()
    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR PEPTIDE DATA, EMITS GFF AND TBI
    //
    PEP_ALIGNMENTS (    reference_tuple,
                        pep_files,
                        max_scaff_size
    )
    ch_versions = ch_versions.mix(PEP_ALIGNMENTS.out.versions)
    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR RNA, NUCLEAR AND COMPLEMENT_DNA DATA, EMITS BIGBED
    //
    GEN_ALIGNMENTS (    reference_tuple,
                        reference_index,
                        gen_files,
                        dot_genome,
                        intron_size
    )
    ch_versions = ch_versions.mix( GEN_ALIGNMENTS.out.versions )
    CDS_ALIGNMENTS (    reference_tuple,
                        reference_index,
                        cds_files,
                        dot_genome,
                        intron_size
    )
    ch_versions = ch_versions.mix( CDS_ALIGNMENTS.out.versions )
    RNA_ALIGNMENTS (    reference_tuple,
                        reference_index,
                        rna_files,
                        dot_genome,
                        intron_size
    )
    ch_versions = ch_versions.mix( RNA_ALIGNMENTS.out.versions )
    emit:
    pep_gff             = PEP_ALIGNMENTS.out.tbi_gff
    gff_file            = PEP_ALIGNMENTS.out.gff_file
    gen_bb_files        = GEN_ALIGNMENTS.out.nuc_alignment
    rna_bb_files        = RNA_ALIGNMENTS.out.nuc_alignment
    cds_bb_files        = CDS_ALIGNMENTS.out.nuc_alignment
    versions            = ch_versions.ifEmpty(null)
```


The .branch() DDL _appears to allow _optional data dependent alignments for peptides, dna, rna and cds inputs._
These in turn involve [PEP_ALIGNMENTS](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/pep_alignments.nf) and
[NUC_ALIGNMENTS ](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/nuc_alignments.nf) subworkflows.
They are complex and need their own decomposition.

They in turn involve miniprot-align [available from the IUC](https://toolshed.g2.bx.psu.edu/repository/browse_repositories?f-free-text-search=miniprot&sort=name&operation=view_or_manage_repository&id=8603bdbca905c70e) in the Toolshed, and PAF2BED and PAFTOOLS that appear to be part of minimap so perhaps supplied from the suite in the Toolshed, plus things already in the Toolshed like samtools_faidx.


```
    include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align/main'
    include { SAMTOOLS_MERGE        } from '../../modules/nf-core/samtools/merge/main'
    include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/samtools/faidx/main'
    include { BEDTOOLS_SORT         } from '../../modules/nf-core/bedtools/sort/main'
    include { BEDTOOLS_BAMTOBED     } from '../../modules/nf-core/bedtools/bamtobed/main'
    include { UCSC_BEDTOBIGBED      } from '../../modules/nf-core/ucsc/bedtobigbed/main'
    include { PAFTOOLS_SAM2PAF      } from '../../modules/nf-core/paftools/sam2paf/main'
    include { PAF2BED               } from '../../modules/local/paf_to_bed'
```


[PAF2BED](https://github.com/sanger-tol/treeval/blob/dev/modules/local/paf_to_bed.nf) is a local module that calls the following command line :


```
paf_to_bed.sh ${file} ${prefix}.bed
```


That bash script is in /treeval/bin so a new tool is needed.

[PAFTOOLS_SAM2PAF](https://github.com/sanger-tol/treeval/blob/dev/modules/nf-core/paftools/sam2paf/main.nf) uses samtools and
perhaps another minimap suite script ? Looks like a new tool is needed:

```samtools view -h ${bam} | paftools.js sam2paf - > ${prefix}.paf
```
