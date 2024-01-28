## [Synteny](https://github.com/sanger-tol/treeval/blob/dev/subworkflows/local/synteny.nf)

This is now implemented as a simple Galaxy subworkflow to make a paf for JBRowse2 viewing.

JBrowse2 synteny tracks have been available since January 20 2024.
Minimap2 will produce a paf.

#### Recipe:

If the reference is `ref` and the syntenic genome reference is `synt`, then 
run minimap2 with `synt` as the reference, and `ref` as the reads to map and be sure to set the advanced options output format to `paf`,
or use the workflow below.

The resulting paf can be provided on the JBrowse tool form together with `synt` as the genome.

Easier to understand as a workflow - the [synteny subworkflow](Galaxy-Workflow-make_synteny_paf_TreeValGal_jan27.ga) will generate a paf:

![image](https://github.com/fubar2/treeval_gal/assets/6016266/29f91b9d-59e8-4a8e-a3d6-b4e9701ef0ff)

When the paf is displayed in JBrowse2, a track like the top panel below will appear - this is the JBrowse peach/grape example
The other panels have a dotplot and a syntenic view, added manually - a few clicks - because cannot figure
out how to automate them :(

![image](https://github.com/fubar2/treeval_gal/assets/6016266/31e8e24a-ea49-44f0-848d-bd296f86d5cf)

For a dotplot track, use the `Add` menu and select `Dotplot view`, as shown below. Make sure the real reference is on the left and the syntenic reference is on the right side. A blank plot usually means the genomes are mixed up for that paf file.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/d0afd4ff-3787-48cb-a542-cd5919fe3bcc)
