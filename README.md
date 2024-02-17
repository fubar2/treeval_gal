# [TreeValGal Project](https://github.com/fubar2/treeval_gal/blob/main/treevalgal/README.md)

*Collaboration between the [Galaxy](https://galaxyproject.org/) and [VGP](https://vertebrategenomesproject.org/) communities to
translate the Sanger [TreeVal NF DDL](https://github.com/sanger-tol/treeval/tree/dev) workflow into something equivalent or better in Galaxy. Everyone with an interest in contributing to
this effort is cordially invited to pitch in.*

## Treevalgal [workflow](treevalgal) 

Depends on JBrowse2 so only available on usegalaxy.eu for testing. 
[Hummingbird sample output](https://usegalaxy.eu/datasets/4838ba20a6d86765a51578280077dbb3/preview) and [Amphyloxus fish sample](https://usegalaxy.eu/datasets/4838ba20a6d8676502c02a6659f467a8/preview) outputs are available. 
This [hummingbird one](https://usegalaxy.eu/datasets/4838ba20a6d86765a82288eace3c126d/preview) has a crude hic track added as a demonstration but it takes a long
time to load a half matrix for a whole chromosome. Probably best turned off unless zoomed in.

These are from a new version that only has a couple of small subworkflows - for making wiggles, hic and paf.

### February 14: For discussion - what's next? 

##### Incorporate tracks from VGP workflow runs into JBrowse2 ?
- What additional tracks would be useful for TreeValGal?
   - for VGP internal use?
   - for scientists and the public?
- What's being generated
  - how to find those tagged VGP WF outputs?
- How to allow access to big JBrowse2 configurations efficiently?
    - JBrowse2 archives contain compressed/indexed reference sequence and track files - so are big.
    - What would these look like to be useful enough to expose with other VGP assembly artifacts?
       - Build centrally with a static public link, for each organism
       - View without redundant data copies using stable URI for remote on-line viewing.

##### Downstream uses for JBrowse2 archives
  - JBrowse2 archive contents can be displayed by a byte-range static web server
    - Typically nginx or apache. 
       - Setup static web pages with links to the unpacked browser archive directory `index.html` file.
  - Can also view on a local laptop browser without Galaxy or internet access
     - a tiny pop-up python webserver is included with the data for local use.
  
##### How best to represent repeat density?
The zoomed in start of this screenshot from the [fish TreeValGal demonstration](https://usegalaxy.eu/datasets/4838ba20a6d8676593a004b88f7c8ab8/preview) shows how a repeat (corresponding to the telomere) counts as 1 for the wiggle on the left, then a group of 3 tiny repeats bump the wiggle to 3. Needs to be normalised to give the proportion of repetitive sequence in a window - the actual count of small repeats might also be useful.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/80331cb3-5459-4c40-b4de-c14a881bf306)

##### TreeValGal modules 

| Module | Status |
|---------------------|-----------------------|
| [yaml_input](yaml_input) | **Not needed** | 
| [gap_finder](gap_finder) | **Prototype available** | 
| [generate_genome](generate_genome) | **Not needed** Existing chromosome lengths tool works in one step. | 
| [hic_mapping](hic_mapping)  | **Fixing hicBuildMatrix for a cool matrix from the VGP hic workflows** | 
| [kmer](kmer)  | **Awaiting fastk and merquryfk tool wrappers**  | 
| [longread_coverage](longread_coverage)  | **Partial prototype available**  | 
| [nuc_alignments](nuc_alignments)  | 
| [pep_alignments](pep_alignments) | 
| [punchlist](punchlist)   | **Need help** - part of hic generation  | 
| [repeat_density](repeat_density)  | **Prototype available.** | 
| [synteny](synteny)  |  **Prototype available.** | 
| [telo_finder](telo_finder) |  **Prototype available in treevalgal workflow now** using seqtk-telo | 

#### January 21

Updating JBrowse1 tool code to JBrowse2 over the past 5 weeks. These are 8 samples - paf and hic are also available.

![image](https://github.com/fubar2/temporary-tools/blob/nohash/jbrowse2/jbrowse8.png)

Working well, but not in IUC yet. Alternative PR is being resurrected so may take some time. 

Deployed on EU as a test tool. See [treevalgal workflow](treevalgal) for more details. Broken wiggle and paf labels in the samples have been fixed.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/7bad8842-3b25-4cfc-a711-29a71eaa6b28)


##### Next steps 
 - Need fastk/merquryFK for a KMER subworkflow.
   - A [Fastk wrapper](https://github.com/galaxyproject/tools-iuc/pull/5550) is currently being developed.
 - Bjoern is thinking about how best to integrate binary HiC format data into the existing infrastructure.
 - Jb2 has specialised Multiple Alignment Format (MAF) tracks, and a blastXML track (translated into GFF3), and those will all be useful in VGP work.
    - awaiting a [PR merge and update](https://github.com/chapmanb/bcbb/pull/141) to BCBIO-GFF conda dependency from Brad Chapman to get these activated.
      
Galaxy's integrated support for genomic feature visualisation at scale will be very hard to match. 


#### December 15
With thanks to Bjoern Gruening and Anna Syme for helping with testing and tools, and the support of Galaxy Australia, there is now a [WIP TreeValGal workflow](treevalgal) on usegalaxy.eu that combines the two gap tracks, the repeats track, a telomere track and the coverage track for the TreeVal small sample test data into a single JBrowse viewer.

Please try it on your own pacbio/refseq data and let me know if this is worth more work to add additional TreeVal tracks to for your use?

Sample images below show how JBrowse does all the work for us. 
All tracks are also in the history as bed files if the user wants them for downstream analyses. 

##### Zoomed out to show windowed bar charts:
![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_gaps_repeats_coverage_zoomout.png)

##### Zoomed midway to show individual features:
![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_repeats_gaps_coverage_midzoom.png)

##### Zoomed in to base level:
![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_gaps_repeats_coverage_basezoom.png)

#### December 13 
1. Some sample outputs from a [prototype longread_coverage workflow](https://github.com/fubar2/treeval_gal/tree/main/longread_coverage#readme) under development and testing. 

#### December 12
1. Some sample outputs from a [prototype repeat_density workflow](https://github.com/fubar2/treeval_gal/tree/main/repeat_density#readme), now available for testing.
Not yet working on the usegalaxy servers. Needs the window_masker tool installed.
   
#### December 10 2023 announcements
1. See sample outputs from a [prototype gap_finder workflow](https://github.com/fubar2/treeval_gal/tree/main/gap_finder#readme), now available for testing.
2. seqtk-telo is now available in the seqtk suite from the Toolshed.

### Work in Progress
Plans are documented for discussion and refinement at [the work plan discussion](https://github.com/fubar2/treeval_gal/discussions/7#discussion-5920941)

For the “rapid” workflow, not all the 18 subworkflows are needed, so this makes a good first target for implementation.
The NF DDL shows the main submodules needed:
```
//
// IMPORT: SUBWORKFLOWS CALLED BY THE MAIN
//
include { YAML_INPUT                                    } from '../subworkflows/local/yaml_input'
include { GENERATE_GENOME                               } from '../subworkflows/local/generate_genome'
include { REPEAT_DENSITY                                } from '../subworkflows/local/repeat_density'
include { GAP_FINDER                                    } from '../subworkflows/local/gap_finder'
include { LONGREAD_COVERAGE                             } from '../subworkflows/local/longread_coverage'
include { TELO_FINDER                                   } from '../subworkflows/local/telo_finder'
include { HIC_MAPPING                                   } from '../subworkflows/local/hic_mapping'
include { KMER                                          } from '../subworkflows/local/kmer'
```
HIC_MAPPING needs some additional subworkflows so the working list for implementation is as above

A first pass at documenting each subworkflow at the level of dependencies and command lines executed at each step has been completed ready for comment and review.
While clarifying the best long term solutions, prototype TreeVal Galaxy subworkflows can be created with existing and some new tools, so users can test them.
These prototypes will help drive the design and features of an optimal production workflow, as discussions progress.

<h2>Background and approach</h2>

Anton has nominated [TreeVal](https://github.com/sanger-tol/treeval/tree/dev) as a test case for translating a NextFlow (NF) workflow into Galaxy.

That translation process needs a clear understanding of the parameters and data TreeVal needs from the user, and what specific application, data and command line is executed at each step.

Once understood, workflows that re-create the same logic can be built using existing or new Galaxy tools, combined into Galaxy subworkflows and a complete Galaxy workflow implemented, tested and documented.

**In the longer term, this is intended to serve as a prototype, to allow users to test it on their own data.**
Reaching consensus on an optimal design when starting from scratch, may be a challenging path, but once all the needed tools are available, Galaxy's GUI workflow editing makes it easy to implement ideas
from users trying things out to incrementally improve working workflows.


Expertise and effort will be needed for:

* NF DDL for insight into how the workflow works
    * Warning! The DDL analysis here is guesswork by Galaxy collaborators, who have never actually used NF!
    * It needs to be reviewed properly.
* the specific scientific context for checking data flows and testing new tools,
* tool and workflow technical skills for building, testing and documentation.

<h3>Sources of truth</h3>

Sanger publish [usage](https://pipelines.tol.sanger.ac.uk/treeval/dev/usage) and [technical ](https://github.com/sanger-tol/treeval/blob/dev/docs/usage.md)
documentation for the [TreeVal](https://github.com/sanger-tol/treeval/tree/dev) workflow source code.
There is an[ overview flow diagram](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_full_diagram.png),
 [detailed output description](https://github.com/sanger-tol/treeval/blob/dev/docs/output.md) for each subworkflow, and an online
 [NF setup generator](https://pipelines.tol.sanger.ac.uk/launch?id=1700725399_4e71a73a94cf).

Test data is available using:

```
curl https://tolit.cog.sanger.ac.uk/test-data/resources/treeval/TreeValTinyData.tar.gz | tar xzf -
```

The tiny test data is 1.2GB


<h3>NF architecture - workflows, subworkflows and modules</h3>

Understanding the NF DDL used to write workflows and modules is key to replicating the functionality needed for Galaxy. Content expertise is also a key ingredient because it can help understand how to adapt NF idioms and idiosyncrasies for Galaxy, confirm process logic, and to test and validate tools and subworkflows as they become available.

NF workflows contain calls to subworkflows and modules. The functionality of a single NF module is usually the Galaxy equivalent of a tool. There are 2 TreeVal workflows . The full one has [18 subworkflows](https://github.com/sanger-tol/treeval/tree/dev/subworkflows/local). They depend on [20 NF-core modules](https://github.com/sanger-tol/treeval/tree/dev/modules/nf-core) and [35 local workflow-specific NF DDL modules](https://github.com/sanger-tol/treeval/tree/dev/subworkflows/local) that depend on [22 specialised scripts and jar files](https://github.com/sanger-tol/treeval/tree/dev/bin) supplied in the _/tree/bin_ directory.  The “rapid” workflow uses only [a subset of those subworkflows](https://github.com/sanger-tol/treeval/blob/dev/workflows/treeval_rapid.nf).

Parameter preparation for module command lines is integrated into the NF workflow DDL. In particular, subworkflows often seem to have specialised DDL involving data and parameters as “channels”, to allow modules to be chained, in contrast to tools connected with noodles in the Galaxy workflow construction GUI.

<h3>Galaxy architecture - workflows, subworkflows and tools</h3>

Galaxy tools are configured to control third-party analysis application code by a wrapper XML document. Applications are hidden from the server behind an abstract tool interface. This isolation layer allows the
Galaxy workflow user to control data and parameter requirements for any tool, using familiar input widgets and explanatory text, without additional programming.
Data flow between any two tools is controlled by connecting an input to an output with matching datatypes in a GUI - no additional programming is required.

As a first try, it makes sense to more or less **re-use the NF subworkflow and module logic in Galaxy**.
Without content expert guidance, there is no information to guide improvements to the existing workflow logic, so the goal is to copy it and duplicate functionality.
Optimisation for Galaxy users then becomes possible because users can respond to the working workflow. It may be possible to build complicated Galaxy tools to perform many steps in a subworkflow, but individual components may be re-used in other parts of related pipelines, so **one tool per module** may be the most adaptable initial approach.

Note that every new tool written will require effort for maintenance. Some will be needed but it is recommended that where it is possible,
the preference will always be to choose an existing, known-good tool for the new subworkflows if it is possible with minor tweaks in the subworkflow itself.

Note: TreeVal uses tabix at the end of many subworkflows. This may be so the outputs can be viewed using Jbrowse. Galaxy's Jbrowse does the compression and indexing automatically so there is no need for that step unless the tabix file is used downstream for some other subworkflow.


<h2>Full TreeVal subworkflows</h2>
The first target is the rapid workflow described above.
Once that is done, new Galaxy subworkflows preserving the logic in
the [full NF Sanger workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_full_diagram.png),
some **_50+ components will be needed_**.
Each more or less corresponds to a Galaxy tool, so the workflow will probably involve 50+ individual tools.

For the full workflow, the subworkflows to be translated into Galaxy are listed below.
As each one is documented, a link will indicate that there is useful material available for that subworkflow.


1. [yaml_input](yaml_input)
2. [ancestral_gene](ancestral_gene)
3. [busco_annotation](busco_annotation)
4. [gap_finder](gap_finder)
5. [gene_alignment](gene_alignment)
6. [generate_genome](generate_genome)
7. [hic_mapping](hic_mapping)
8. insilico_digest
9. [kmer](kmer)
10. [longread_coverage](longread_coverage)
11. [nuc_alignments](nuc_alignments)
12. [pep_alignments](pep_alignments)
13. pretext_ingestion
14. [punchlist](punchlist)
15. [repeat_density](repeat_density)
16. selfcomp
17. [synteny](synteny)
18. [telo_finder](telo_finder)
    


