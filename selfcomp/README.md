# treeval_gal
Sanger treeval nf workflow translation into Galaxy work in progress

![Flow chart](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_self_comp.png)


The selfcomp subworkflow is a comparative genomics analysis algorithm originally performed by the Ensembl projects database, and reverse engineered in Python3 by @yumisims. It involves comparing the genes and genomic sequences within a single species. The goal of the analysis is to identify haplotypic duplications in a particular genome assembly.
Output files

    treeval_upload/
        *_selfcomp.bigBed: BigBed file containing selfcomp track data.

