## TreeValGal demonstration workflow

This workflow integrates tracks from prototype and in development TreeVal subworkflows into a single JBrowse configuration ready to view. 
Testing results, suggestions and contributions are very welcome.

![image](https://github.com/fubar2/treeval_gal/assets/6016266/7c74faef-3cf3-490e-8ebe-fcd55687e99c)

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
