## TreeValGal workflow 

### May 2: Core TreeValGal components

See the main [README.md](README.md) for the latest updates to the workflow and subworkflows.

The TreeValGal workflow is complex and overwhelming at first glance, so the major components, including 3 subworkflows, are described here.

The JBrowse2 tool on the right side is the most prominent workflow element, producing an integrated interactive JBrowse2 browser configured with 15 tracks. 
Some tracks (the paf, gaps and telomeres) are turned off in the default view to reduce the complexity at opening the browser screen. 
All tracks can be turned on and off using the JBrowse2 track menu to allow the user to focus on tracks of specific interest without distraction.

#### 1. Workflows components reporting the fasta reference.

These do not involve external data other than PacBio reads fastqsanger input, and they produce four JBrowse2 tracks for display:
 - Repeatmasker repeats as a GFF3 track
 - Windowmasker repeats as a bigwig track
 - Gaps as a bigwig track
 - PacBio depth of coverage as a bigwig track

![image](https://github.com/fubar2/treeval_gal/assets/6016266/5fc8697a-ae5a-42f2-9536-99b1a35aa309)

The top row of tools and data flows run the tool `gfastats` on the reference fasta to report sequence gaps as a bigwig track. The second row 
prepares a file with window start/end coordinates, used to count features as a bigwig track for gaps, coverage and repeats. The third row runs the model-free `Windowmasker` 
repeat finder tool and the fourth uses `Minimap2` to map PacBio reads to the reference for coverage. The `bamcoverage` tool creates a bigwig directly from the mapped bam file.

The `repeatmasker` tool seen at the bottom of the workflow section shown above, creates a `GFF3` track of repeats. A DFam taxonId can be provided to over-ride the default (*Homo sapiens*) taxon.
The masked fasta output is retained for use in the creation of sequence similarity paf files from `mashmap` 
to minimise noise from uninformative low-complexity sequence matches.

Two bigwig tracks (gaps and windowmasker repeats) use a specialised subworkflow to count depth of coverage over the generated window start/end coordinates based on the user window size (default is 100bp) input. 

![TreeValGal bed to bigwig subworkflow](https://github.com/fubar2/treeval_gal/assets/6016266/146176b3-27c0-4b82-8bac-11d1959933ce)

#### 2. Workflow components reporting external fasta annotation data.

The middle segment of the main TreeValGal workflow canvas is responsible for mapping NCBI or other external annotation fasta files to create JBrowse2 `GFF3` and `bed` tracks.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/e9a88a46-5901-4851-8ed3-436d8683c75b)

The subworkflow for these tracks is described in [detail here](../gene_alignment/README.md)


#### 3. Workflow components preparing pairwise mapping format (PAF) tracks for haplotype self-comparison and closely related species sequence synteny or similarity.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/db642299-eaa8-44fc-986a-243841809214)

This process uses the `mashmap` tool multiple times with a range of settings, to provide choices for the user, depending on how closely related the provided genomes are, in a third subworkflow.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/1ca4da9f-b0b4-4e85-8c23-656932185e14)

Like the [gene_alignment equivalent](../gene_alignment/README.md), this workflow uses `pick` workflow logic components that make all the `mashmap` steps optional, since the `mashmap` tool
does not natively allow optional inputs.


### April 23 update

The latest TreeValGal workflow is on the EU server at [https://usegalaxy.eu/published/workflow?id=2e93613c688ea52e](https://usegalaxy.eu/published/workflow?id=2e93613c688ea52e)


### February 21 update

The [TreeValGal WF](Galaxy-Workflow-TreeValGal_Feb11.ga) is [available on usegalaxy.eu](https://usegalaxy.eu/published/workflow?id=4d44b055b557cba5) depends on JBrowse2 so only available on usegalaxy.eu for testing.

[Hummingbird sample output](https://usegalaxy.eu/datasets/4838ba20a6d8676565fed5852c79ff4d/preview) and [Amphioxus fish sample](https://usegalaxy.eu/datasets/4838ba20a6d867655aedf35d84ed3d59/preview) outputs are available. 

These JBrowse configurations now include repeatmasker GFF tracks, from the latest Feb_11 revision, that only has a couple of small subworkflows - for making wiggles and for optional hic and paf.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/aa6c2e1b-e7c8-4149-8746-5c710b58afb4)

The wiggle maker is the most complicated subworkflow and is used for 3 tracks.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/75bdf6aa-62af-4cbf-a8e2-d650ae4314d1)

The optional synteny and hic track subworkflows are relatively trivial...


### January 2024 update

This workflow integrates tracks from prototype TreeVal subworkflows into a single JBrowse2 configuration, ready to view, share and download.

> [!WARNING]
> The compressed archives can be very big, because they contain the reference sequences and tracks, albeit compressed and indexed.

Testing results, suggestions and contributions are very welcome.

Using the hummingbird test data from Anna Syme on the EU server, it produces [this live JBrowse2 instance](https://usegalaxy.eu/datasets/4838ba20a6d86765abb0a6a22fe9087d/preview). Nadolina's Lancet fish shown below was chosen as the synteny genome. 

Only a solitary telomere. Does the bird need a different telomere repeat?

Fixes for the uninformative wiggle and paf track names will be live on EU shortly....

Sample image after manually adding a dot plot:
![image](https://github.com/fubar2/treeval_gal/assets/6016266/8724c707-8fa3-4b34-8b04-708d93c9a28e)


Clicking on a syntenic feature shows the details of the match and the fish sequence if wanted:

![image](https://github.com/fubar2/treeval_gal/assets/6016266/e4bbf907-8555-45d8-94ba-ddb57040760e)

Nadolina's Lancet fish - same workflow:

![image](https://github.com/fubar2/treeval_gal/assets/6016266/0ace5ce1-e5d4-4f34-864a-0f0bbaa27bb5)

and this time, the detailed view of a syntenic region shows a syntenic segment of bird sequence
![image](https://github.com/fubar2/treeval_gal/assets/6016266/dace8f97-ee02-4eed-bd85-7a615a5849f3)



The main treevalgal workflow and the subworkflows it calls are shared on EU as `treevalgal_jan27` from **fubar** :

![image](https://github.com/fubar2/treeval_gal/assets/6016266/2811a123-d073-4128-b96b-058ca72c79ae)


#### December 15 prototype

With thanks to Bjoern Gruening and Anna Syme for help with testing and tools, and the support of Galaxy Australia, 
the current (December 15) version combines the two gap tracks, the repeats track and the coverage track for the TreeVal small sample test data.

*Note that there are only 3 gaps in the entire test reference so hard to find them, and the pacbio sample only covers some of a single contig, so select ENA|OV656687|OV656687.1 for display, otherwise there will not be any coverage to see in this tiny sample demonstration.*

A demonsration history with the viewable preconfigured JBrowse is shared here [on usegalaxy.eu](https://usegalaxy.eu/u/fubar/h/vgpdemogapsrepeatscoveragecombineddec15) and so is [the prototype workflow](https://usegalaxy.eu/u/fubar/w/vgpdemogapsrepeatscoveragecombined) so please try it on your own pacbio/refseq data.

Sample images show how JBrowse does all the work of density and other displays based on the zoom level. 
All tracks are also in the history as bed files if the user wants them for downstream analyses. 

##### Zoomed out to show windowed bar charts:

![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_gaps_repeats_coverage_zoomout.png)

##### Zoomed midway to show individual features:

![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_repeats_gaps_coverage_midzoom.png)

##### Zoomed in to base level:

![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_gaps_repeats_coverage_basezoom.png)
