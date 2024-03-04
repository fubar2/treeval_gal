# [TreeValGal Project](https://github.com/fubar2/treeval_gal/blob/main/treevalgal/README.md)

*Collaboration between the [Galaxy](https://galaxyproject.org/) and [VGP](https://vertebrategenomesproject.org/) communities to
translate the Sanger [TreeVal NF DDL](https://github.com/sanger-tol/treeval/tree/dev) workflow into something equivalent or better in Galaxy. Everyone with an interest in contributing to
this effort is cordially invited to pitch in.*

Origins and early progress [documented back here](#December15)

### Treevalgal workflow output samples

[Hummingbird sample output](https://usegalaxy.eu/datasets/4838ba20a6d86765cecca865d35c5146/preview) and [Amphioxus fish sample](https://usegalaxy.eu/datasets/4838ba20a6d86765369b13bd5c7199ce/preview) outputs are available. 
This [hummingbird one](https://usegalaxy.eu/datasets/4838ba20a6d86765a82288eace3c126d/preview) has a crude hic track added as a demonstration but it takes a long
time to load a half matrix for a whole chromosome. Probably best turned off unless zoomed in. Not needed since the current pretextview solution works well.

These are from [the current version](treevalgal) that only has a couple of small subworkflows - for making wiggles, hic and paf.

TreeValGal depends on the updated [Jbrowse2 tool](#JBrowse2) to present repeats, gaps and other features as tracks in a configurable browser. 
JBrowse2 is currently only available on `usegalaxy.eu` for testing. 

### 4 March: Updated work plans

- Anna working on nucleotide and peptide alignments subworkflows
  - ncbi downloader improvements suggested
- Next JB2 release has optional inputs for workflows. 
- New AutogenJB2 (collection -> JB2) tool soon available on EU
  - The core of the JB2 track generator requires only a path, name and ext in most cases
    - paf will need a subcollection of paf, reference1...referenceN

### 25 Feb: Big picture ideas, mashmap and next steps

- Long term goal
  - Add interactive JBrowse configurations to genomeArk web pages, to make VGP data easily accessible as tracks to anyone with a web browser.
  
- Intermediate milestones
   - Production TreeValGal for internal VGP use
   - Build and distribute other VGP JB2 track sets for different audiences
   - Provide a test bed for the development and display of automated annotation
     - integrate useful tracks into TreeValGal and other outputs as they become available.
   - Additional JBrowse2 GTN training modules
       - create dot plots and reverse mappings from mashmap PAF tracks ... 
    
- Short term milestones
  - a "good to test" track list with tool tuning parameters as needed.
  - Implement as the initial TreeValGal workflow and hand over
  - Clean-up and acceptance-testing as part of a routine VGP WF run. 
  - Busco genes could be added to TreeValGal as GFF3 if available
    - Action: could someone please make a sample busco gff3 for one of the test assemblies we are currently using: calanna and bralanc ?
    - a gff is probably more suited for a JB2 track than a bed
      - easy enough to include both if useful  
      
- The next release (_5) of JBrowse2 on EU will optionally create a zip file instead of being a viewable HTML page. 
  - This enables the Galaxy uploader to push a JBrowse directory to S3. It seemed confused by html
    - Public website proof of concept using a zip from Galaxy to an S3 bucket thanks to Bjoern.
      
- Nadolina has nominated the [gene_alignment](gene_alignment) TreeVal subworkflow to implement.
    - Tools should all be on EU but may have missed something.
    - Two main subworkflows are used - for nucleotide and peptide mapping so each of those can be built and tested
    - Once tested, we can figure out how best to fit in with the existing WF.
   
- mashmap proposed as a way to make multi-organism PAF tracks.
  - Ross built and tested a prototype
  - Bjoern has take the tool over, and made it far better.
  - Will soon be available from IUC
  - Sample 3 way PAF display below was manually set up, starting with a JBrowse with a single mashmap track and arabidopsis reference.
    - Other two genomes added as separate linear tracks below the main reference gives reverse mappings!
    - Easy to configure with JB2 - a few moments clicking - all the references are available.
    - Dot plots too - all there but require manual exposure
    
![image](https://github.com/fubar2/treeval_gal/assets/6016266/fae539ba-55e0-4fa3-b56c-3e5b96027e3f)


Here are some bird on bird dot plots - they even look a bit like a flock of birds if you squint.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/34e10961-4d55-4bc0-8342-295e03333976)

### 21 Feb: For discussion - if it looks potentially useful, what's next?

#### 1. Add tracks and other organisms ?
  - What additional tracks would be useful?
     - for VGP internal use? 
       - Other TreeVal modules: selfcomp
          - kmer/busco/pretext etc already in VGP)
     - for external researchers
     - for the public - making VGP data really easy to "access"
  - Track files being generated routinely can be added if there's a URI or path.
    - Which of [the existing VGP WF](https://galaxyproject.org/projects/vgp/workflows/) outputs would be useful to add?
      - tracks take some disk space but there's no real technical limit to the size of a track menu. 
  - Other VGP assemblies to run TreeValGal on to serve as samples?
    
#### 2. Jbrowse2 outputs are self contained and portable
  - JBrowse2 directories support interactive display
    -  Any byte-range request static web server like nginx or apache. 
  - Can also view on a local laptop browser without Galaxy or internet access
     - a tiny pop-up python webserver `jb2_webserver.py` is included with the data for local displays.
  - How to allow access to big JBrowse2 configurations efficiently?
    - JBrowse2 archives contain compressed/indexed reference sequence and track files - so are big.
  
#### 3. How best to represent repeat density?
Anna asked about a repeatmasker version of the repeats, so there's a gff track (no dfam species and defaults) to compare to the windowmasker bed and wiggle on the current [Amphioxus JBrowse2](https://usegalaxy.eu/datasets/4838ba20a6d867655aedf35d84ed3d59/preview).

![image](https://github.com/fubar2/treeval_gal/assets/6016266/860af67f-00bc-4848-af2f-ba47e9762710)

Repeatmasker seems more selective with half or fewer but bigger hits. It is model based so takes much longer (~2 orders of magnitude) to complete although that's not important for a workflow. 7 hours 4 cores vs 15 minutes 1 core for the hummingbird for example.

The other repeat data is from a model-free method - [NCBI windowmasker](https://doi.org/10.1093/bioinformatics/bti774) - so no surprise that a model-based method like repeatmasker is verrrrry different with defaults - suggestions appreciated!

Here's a zoomed in part of the current hummingbird TreeValGal output after hiding the other tracks and opening the repeatmasker gff and the repeat bed. The wiggle is windowmaker smoothed over 1k windows and reflects that detail differently, particularly when zoomed in.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/3af7e2a8-e2c7-4e32-995d-6619435caaf9)

Each method and view gives different information about the complex underlying biology here and JBrowse2 makes it easy to pick and choose between them.

#### 4. Current track list
![image](https://github.com/fubar2/treeval_gal/assets/6016266/8daa9117-9a5e-4637-bb62-8bd921809813)


### Potential new resource?

Generating and publishing updated live JBrowse2 links for any VGP page is possible, because JBrowse2 configuration can be completely machine driven. 
- Offer convenient, highly granular tracks derived from VGP and other data
- readily accessible to anybody with internet, a web browser and a question. 
- Adding automated feature annotation tracks could build a potentially valuable, accessible and sustainable resource for genome sciences.

Each JBrowse2 directory is independent, with an `index.html`, track files and javascript.
- served typically by Apache or Nginx
- potentially stored (untested!) on an S3 bucket.

Most of the components are available. 
- An initial functional track menu can quickly be made available.
- It will require storage and effort TBA
- Proof of concept sample outputs below. 

Challenges.
- storage for compressed/indexed data (a few GB each perhaps)
- computational resources required are very substantial for creation.
- Individual updates are only needed when data changes.
- Workflow change means recreating everything
 - not a problem if "reuse jobs" really works...
- ongoing dedicated effort for sustainable development.
  
Given those, configurations can be improved, by editing the workflows, to suit real user needs.
The [rat database does this](https://rgd.mcw.edu/wg/help3/tools/rgd-genome-browsers/the-rat-jbrowse-genome-browser/) with JB1, so it's not novel, but would add a lot of value to existing data with relatively little additional resources.


### JBrowse2

![image](https://github.com/fubar2/treeval_gal/assets/6016266/8e03a285-efba-49e9-840e-cf9dbfb6b4ec)

Now does 10 track types - any of these datatypes in a Galaxy history or workflow can be added as a track to a JBrowse2 instance. Outputs are immediately viewable in Galaxy. 
The JBrowse2 directory can be downloaded as a zip file, for local unzipping and viewing using a provided tiny python web server. 
These can also be served using any byte-range request static web server. 
They tend to be huge - GB in size, with the reference and all tracks in compressed and indexed formats.

    bam
    bed
    bigwig
    blastxml to gff3
    cram
    gff3
    hic
    maf
    paf
    vcf

### TreeValGal modules 

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

### January 21

Updating JBrowse1 tool code to JBrowse2 over the past 5 weeks. These are 8 samples - paf and hic are also available.

![image](https://github.com/fubar2/temporary-tools/blob/nohash/jbrowse2/jbrowse8.png)

Working well, but not in IUC yet. Alternative PR is being resurrected so may take some time. 

Deployed on EU as a test tool. See [treevalgal workflow](treevalgal) for more details. Broken wiggle and paf labels in the samples have been fixed.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/7bad8842-3b25-4cfc-a711-29a71eaa6b28)


#### Next steps 
 - Need fastk/merquryFK for a KMER subworkflow.
   - A [Fastk wrapper](https://github.com/galaxyproject/tools-iuc/pull/5550) is currently being developed.
 - Bjoern is thinking about how best to integrate binary HiC format data into the existing infrastructure.
 - Jb2 has specialised Multiple Alignment Format (MAF) tracks, and a blastXML track (translated into GFF3), and those will all be useful in VGP work.
    - awaiting a [PR merge and update](https://github.com/chapmanb/bcbb/pull/141) to BCBIO-GFF conda dependency from Brad Chapman to get these activated.
      
Galaxy's integrated support for genomic feature visualisation at scale will be very hard to match. 


### December15
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
    


