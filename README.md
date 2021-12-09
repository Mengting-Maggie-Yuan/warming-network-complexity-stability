This repository contains code used to generate results and figures for the manuscript "Climate Warming Enhances Microbial Network Complexity and Stability".
Citation: Yuan, M.M., Guo, X., Wu, L. et al. Climate warming enhances microbial network complexity and stability. Nat. Clim. Chang. (2021). https://doi.org/10.1038/s41558-021-00989-9

Part of the analyses in this paper is available through the <b>Molecular Ecological Network Analysis Pipeline (MENAP)</b>: http://ieg4.rccc.ou.edu/MENA/. Results generated in MENAP is noted in the following instruction.

## Updates

2021.12.8 Bug fix in FigS1.LTED/FigS1ab.find_env_links_example.sh and FigS1.LTED/FigS1ab.find_env_links.py  


## Instructions for generating result/figures in this manuscript

- Figure 1: MENAP - Global Network properties, Randomize the network structure and then calculate network properties, or through R package igraph
- Figure 2: MENAP - Global Network properties, Randomize the network structure and then calculate network properties, or through R package igraph
- Figure 3: See folder Fig3_and_S7.stability
- Figure 4: Available in R package vegan

- Extended Data Figure 1: Available in R package vegan and agricolae
- Extended Data Figure 2: Available through original study introduced cohesion - Herren, C., McMahon, K. Cohesion: a method for quantifying the connectivity of microbial communities. ISME J 11, 2426–2438 (2017). https://doi.org/10.1038/ismej.2017.91
- Extended Data Figure 3: See folder ExtendedDataFig3.module_preservation
- Extended Data Figure 4: MENAP - Output for the Cytoscape software visualization, or through R package igraph
- Extended Data Figure 5: See folder ExtendedDataFig5.constancy

- Supplementary Figure 1: See folder FigS1.LTED
- Supplementary Figure 2: Available through original study - Goberna, M, Montesinos‐Navarro, A, Valiente‐Banuet, A, et al. Incorporating phylogenetic metrics to microbial co‐occurrence networks based on amplicon sequences to discern community assembly processes. Mol Ecol Resour. 2019; 19: 1552– 1564. https://doi.org/10.1111/1755-0998.13079
- Supplementary Figure 3: Available through R
- Supplementary Figure 4: Differential abundance analysis is available through the original study - Morton, J.T., Marotz, C., Washburne, A. et al. Establishing microbial composition measurement standards with reference frames. Nat Commun 10, 2719 (2019). https://doi.org/10.1038/s41467-019-10656-5, and Songbird software - https://github.com/biocore/songbird
- Supplementary Figure 5: Available through R package igraph
- Supplementary Figure 6: MENAP - Module-EigenGene analyses
- Supplementary Figure 7: See folder Fig3_and_S7.stability
- Supplementary Figure 8: Available in R package vegan
- Supplementary Figure 9: MENAP - Global Network properties, or through R package igraph
- Supplementary Figure 10: FastSpar - Stephen C Watts, Scott C Ritchie, Michael Inouye, Kathryn E Holt, FastSpar: rapid and scalable correlation estimation for compositional data, Bioinformatics, Volume 35, Issue 6, 15 March 2019, Pages 1064–1066, https://doi.org/10.1093/bioinformatics/bty734, and MENAP - Global Network properties, or through R package igraph

## Instructions for input data used in the code

- Fig3_and_S7.stability  
<b>CorrelationMatrix_Y14_W.txt</b>: correlation matrix downloaded from MENAP, an upper triangle showing correlation values of OTU pairs, without header row or column names. The order of OTUs can also be downloaded from MENAP. Note that only OTUs occur more than 12 in 24 samples are included in correlation calculation in this study. Sample data is from Year 2014 warming plots.  
<b>OTUtable_NetworkedOTUs_AllSamples.txt</b>: An OTU table containing the abundances of OTUs only occurred in each network. If an OTU is absent from a network, its abundance is set to zero in all samples from that year/treatment.  
<b>SampleMap_AllSamples.txt</b>: A metadata file containing treatment information of the samples in OTUtable_NetworkedOTUs_AllSamples.txt. Sample order in these two files should match.  
<b>OTUtable_AllOTUs_Y14_W.txt</b>: The original OTU table containing all the detected OTUs before construction of the network. Sample data is from Year 2014 warming plots.  
<b>NodeAttribute_Y14_W.txt</b>: Network node attribute table. Can be downloaded through MENAP. Sample data is from Year 2014 warming plots.  

- ExtendedDataFig3.module_preservation  
<b>NetworkNodeTable_AllNetworks.txt</b>: Network node attribute table from all the networks.  

- ExtendedDataFig5.constancy  
<b>OTUtable_NetworkedOTUs_AllSamples.txt</b>: An OTU table containing the abundances of OTUs only occurred in each network. If an OTU is absent from a network, its abundance is set to zero in all samples from that year/treatment.  
<b>SampleMap_AllSamples.txt</b>: A metadata file containing treatment information of the samples in OTUtable_NetworkedOTUs_AllSamples.txt. Sample order in these two files should match.  
<b>NetworkNodeTable_AllNetworks.txt</b>: Network node attribute table from all the networks.  
<b>NetworkEdgeTable_AllNetworks.txt</b>: Network edge attribute table from all the networks.  

- FigS1.LTED  
<b>OTUtable_AllOTUs_Y14_W.txt</b>: The original OTU table containing all the detected OTUs before construction of the network. Sample data is from Year 2014 warming plots.  
<b>EnvironemtnalVariables_Y14_W.tsv</b>: Environmental variable values for each sample. Sample data is from Year 2014 warming plots.  
<b>Distance_Y14_W.tsv</b>: Two dimensional geodistance in meters, from the center of each plot to a defined origin (in this case the center of Plot 1 in our field experiment setting. Sample data is from Year 2014 warming plots.  
<b>NetworkEdgeTable_Y14_W.txt</b>: Network edge attribute table. Can be downloaded through MENAP. Sample data is from Year 2014 warming plots.  
