## [Synteny](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/synteny.nf)

This is easy as a simple subworkflow to make a paf for viewing.

JBrowse2 synteny tracks have been available since January 20 2024.
Minimap2 will produce a paf.

#### Recipe:

If the reference is ref and the syntenic genome reference is synt, then 
run minimap2 with `synt` as the reference, and `ref` as the reads to map. Set the output to `paf`

The resulting paf can be provided on the JBrowse tool form together with `synt` as the genome.

A track like the top panel below will appear - this is the JBrowse peach/grape example
The other panels have a dotplot and a syntenic view, added manually - a few clicks - because cannot figure
out how to automate them :(

![image](https://github.com/fubar2/treeval_gal/assets/6016266/31e8e24a-ea49-44f0-848d-bd296f86d5cf)



