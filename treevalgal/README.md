## TreeValGal demonstration workflow January 2024 update

This workflow integrates tracks from prototype TreeVal subworkflows into a single JBrowse2 configuration, ready to view, share and download.

> [!WARNING]
> The compressed archives can be very big, because they contain the reference sequences and tracks, albeit compressed and indexed.

Testing results, suggestions and contributions are very welcome.

Using the hummingbird test data from Anna Syme on the EU server, it produces [this live JBrowse2 instance](https://usegalaxy.eu/datasets/4838ba20a6d86765a3a465ade83523e0/preview). Nadolina's Lancet fish shown below was chosen as the synteny genome. 

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
