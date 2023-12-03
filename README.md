# treeval_gal

### tl;dr

Translate the Sanger [TreeVal NF DDL](https://github.com/sanger-tol/treeval/tree/dev) workflow into Galaxy. Extends ongoing collaboration between the [Galaxy](https://galaxyproject.org/) and [VGP](https://vertebrategenomesproject.org/) communities. Just started. If you would like to contribute, please join us. 

### Work in Progress - documenting TreeVal

For the “rapid” workflow, not all the 18 subworkflows are needed, so this makes a good first target for implementation.
Links indicate that the subworkflow logic and command lines have been documented. There may be others that are called by ones listed here - will
only know when all these have been documented.

1. [yaml_input](yaml_input)
4. [gap_finder](gap_finder)
6.  [generate_genome](generate_genome)
7. [hic_mapping](hic_mapping)
9. [kmer](kmer)
10. longread_coverage
11. [nuc_alignments](nuc_alignments)
12. pep_alignments
15. [repeat_density](repeat_density)
18. [telo_finder](telo_finder)


<h2>Background</h2>

Anton has nominated [TreeVal](https://github.com/sanger-tol/treeval/tree/dev) as a test case for translating a NextFlow (NF) workflow into Galaxy.

The first challenge is to understand and document what parameters and data TreeVal needs from the user, and what it does at each step.

Steps can be re-created using existing or new Galaxy tools, combined into Galaxy subworkflows and a complete Galaxy workflow implemented, tested and documented.

Dedicated resources including expertise and effort will be needed for:



* NF DDL for insight into how the workflow works,
* the specific scientific context for checking data flows and testing new tools,
* tool and workflow technical skills for building, testing and documentation.

<h3>Resources</h3>


Sanger publish [usage](https://pipelines.tol.sanger.ac.uk/treeval/dev/usage) and [technical ](https://github.com/sanger-tol/treeval/blob/dev/docs/usage.md)documentation for the [TreeVal](https://github.com/sanger-tol/treeval/tree/dev) workflow source code. There is an[ overview flow diagram](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_full_diagram.png), [detailed output description](https://github.com/sanger-tol/treeval/blob/dev/docs/output.md) for each subworkflow, and an online [analysis service](https://pipelines.tol.sanger.ac.uk/launch?id=1700725399_4e71a73a94cf).

Test data is available using:


```
curl https://tolit.cog.sanger.ac.uk/test-data/resources/treeval/TreeValTinyData.tar.gz | tar xzf -
```


The download is 1.2GB so the tiny test data is vast.

<h3>NF workflow DDL code</h3>


Understanding the NF DDL used to write workflows and modules is key to replicating the functionality needed for Galaxy. Content expertise is also a key ingredient because it can help understand how to adapt NF idioms and idiosyncrasies for Galaxy, confirm process logic, and to test and validate tools and subworkflows as they become available.

There will be new tools, mostly trivial bash or other “one liners” but some more complex. Once the inputs, parameters and processing are well understood, the components available as Galaxy tools and the subworkflows as Galaxy subworkflows, then it will be possible to put the entire thing together for testing and documentation.

<h3>NF architecture - workflows, subworkflows and modules</h3>


NF workflows contain calls to subworkflows and modules. The functionality of a single NF module is usually the Galaxy equivalent of a tool. There are 2 TreeVal workflows . The full one has [18 subworkflows](https://github.com/sanger-tol/treeval/tree/dev/subworkflows/local). They depend on [20 NF-core modules](https://github.com/sanger-tol/treeval/tree/dev/modules/nf-core) and [35 local workflow-specific NF DDL modules](https://github.com/sanger-tol/treeval/tree/dev/subworkflows/local) that depend on [22 specialised scripts and jar files](https://github.com/sanger-tol/treeval/tree/dev/bin) supplied in the _/tree/bin_ directory.  The “rapid” workflow uses only [a subset of those subworkflows](https://github.com/sanger-tol/treeval/blob/dev/workflows/treeval_rapid.nf).

Parameter preparation for module command lines is integrated into the NF workflow DDL. In particular, subworkflows often seem to have specialised DDL involving data and parameters as “channels”, to allow modules to be chained, in contrast to tools connected with noodles in the Galaxy workflow construction GUI.

<h3>Galaxy architecture - workflows, subworkflows and tools</h3>


In Galaxy workflows, tools control third-party analysis application code through an abstract interface, configured by the wrapper XML document. This isolation layer allows the Galaxy server to automatically provide data and parameter requirements in a GUI, using familiar input widgets and explanatory text. It also allows Galaxy tools to interoperate transparently for users as they connect them with noodles in workflows.

As a first try, it makes sense to more or less **re-use the NF subworkflow and module logic in Galaxy**.
Without content expert guidance, there is no information to guide improvements to the existing workflow logic, so the goal is to copy it and duplicate functionality.
Optimisation for Galaxy users then becomes possible because users can respond to the working workflow. It may be possible to build complicated Galaxy tools to perform many steps in a subworkflow, but individual components may be re-used in other parts of related pipelines, so **one tool per module** may be the most adaptable initial approach.

Note that every new tool written will require effort for maintenance. Some will be needed but it is recommended that where it is possible,
the preference will always be to choose an existing, known-good tool for the new subworkflows if it is possible with minor tweaks in the subworkflow itself.


<h2>NF component descriptions and work plans</h2>

For new Galaxy subworkflows preserving the logic in
the [full NF Sanger workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/v1-1-0/treeval_1_1_0_full_diagram.png),
some **_50+ components will be needed_**.
Each more or less corresponds to a Galaxy tool, so the workflow will probably involve 50+ individual tools.
Some are already available in the Toolshed but many are not.

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
10. longread_coverage
11. nuc_alignments
12. pep_alignments
13. pretext_ingestion
14. punchlist
15. [repeat_density](repeat_density)
16. selfcomp
17. synteny
18. [telo_finder](telo_finder)

The DDL for these subworkflows must be decomposed into their steps, data flows, and transformation code resources. Module steps must be implemented as Galaxy tools and validated by a content expert.

That is what this repository is for.
