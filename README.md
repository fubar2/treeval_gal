# [treevalgal](https://github.com/fubar2/treeval_gal/blob/main/treevalgal/README.md)

*New proposal in the ongoing collaboration between the [Galaxy](https://galaxyproject.org/) and [VGP](https://vertebrategenomesproject.org/) communities to
translate the Sanger [TreeVal NF DDL](https://github.com/sanger-tol/treeval/tree/dev) workflow into Galaxy. Everyone with an interest in contributing to
this effort is cordially invited to pitch in.*

## Current status of TreeVal modules needed for the rapid workflow 

1. [yaml_input](yaml_input) **Not needed**
4. [gap_finder](gap_finder) **Prototype available**
6. [generate_genome](generate_genome) **Not needed** Existing chromosome lengths tool works in one step.
7. [hic_mapping](hic_mapping)
9. [kmer](kmer)
10. [longread_coverage](longread_coverage) **Partial prototype available**
11. [nuc_alignments](nuc_alignments)
12. [pep_alignments](pep_alignments)
14. [punchlist](punchlist) **Redundant** - avoiding PAF removes the need for complication
15. [repeat_density](repeat_density) **Prototype available.**
18. [telo_finder](telo_finder) **Prototype available in treevalgal workflow now** using seqtk-telo

## News

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
17. synteny
18. [telo_finder](telo_finder)

