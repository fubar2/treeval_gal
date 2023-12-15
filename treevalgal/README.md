## TreeValGal demonstration workflow

![image](https://github.com/fubar2/treeval_gal/assets/6016266/7c74faef-3cf3-490e-8ebe-fcd55687e99c)

With thanks to @bgruening and @annasyme for helping with testing and tools, there is now a WIP workflow on usegalaxy.eu that combines the two gap tracks, the repeats track and the coverage track for the TreeVal small sample test data into a single JBrowse viewer.

*Note that there are only 3 gaps in the entire test reference so hard to find them, and the pacbio sample only covers some of a single contig, so make sure you select ENA|OV656687|OV656687.1 for display, otherwise there is no coverage to show in this demonstration.*

[History is shared on usegalaxy.eu](https://usegalaxy.eu/u/fubar/h/vgpdemogapsrepeatscoveragecombineddec15) and so is [the prototype workflow](https://usegalaxy.eu/u/fubar/w/vgpdemogapsrepeatscoveragecombined) so please try it on your own pacbio/refseq data.

Goal is to provide a demonstration by adding tracks from other TreeVal subworkflows for your use, 
Suggestions and contributions are very welcome.

Sample images show how JBrowse does all the work for us. 
All tracks are also in the history as bed files if the user wants them for downstream analyses. 

##### Zoomed out to show windowed bar charts:
![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_gaps_repeats_coverage_zoomout.png)

##### Zoomed midway to show individual features:
![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_repeats_gaps_coverage_midzoom.png)

##### Zoomed in to base level:
![image](https://github.com/fubar2/treeval_gal/blob/main/vgp_gaps_repeats_coverage_basezoom.png)
