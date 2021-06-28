# Low-Pass Whole Genome Bisulfite Sequencing of Neonatal Dried Blood Spots

[![DOI](https://zenodo.org/badge/291577655.svg)](https://zenodo.org/badge/latestdoi/291577655)

Each folder in this repository contains the scripts required to recreate the analyses and plots for the main figures from the manuscript: "Low-Pass Whole Genome Bisulfite Sequencing of Neonatal Dried Blood Spots Identifies a Role for RUNX1 in Down Syndrome DNA Methylation Profiles". 

## Scripts

- `DMR_node.sh`: Example of how to perform multiple DM.R runs in parallel on a login node. 

- `Figure 1`: Distinct DS methylome profiles.
	- `cnvPlot.R`: Plot of chromosome 21 copy number using output from [CNV_Me](https://github.com/hyeyeon-hwang/CNV_Me).
	- `global.R`: Process all samples for equal CpG coverage and smooth, test for global methylation differences, and create a single CpG density plot.
	- `heatmapUpdate.R`: Update DMRichR heatmaps from multiple comparisons to have a consistent style.
	
- `Figure 2`: Gene Ontology enrichments.
	- `GOplot.R`: Slim and compare GOfuncR terms from one comparison to other comparisons.
	
- `Figure 3`: Consensus DMR profiles.
	- `consensusDMRs.sh`: Create a consensus DMR profile from multiple comparisons for `consensusPCA.R`.
	- `ConsensusPCA.R`: PCA of all samples using consensus DMRs.
	- `machineLearning.R`: DMR machine learning feature selection analysis and heatmap.
	- `venn.R`: Create an area-proportional Venn diagram based on overlaps of genomic coordinates from multiple comparisons.
	
- `Figure 4`: Divergent DNA hyper- and hypo-methylation profiles.
	- `chromatinPlots`: Parsing and plots for Roadmap epigenomics chromHMM core 15-state model and 5 core histone modifications from 127 reference epigenomes.
		- `chromatinPlots.R`: Parsing and summary heatmaps for multiple LOLA chromHMM and core histone modification enrichment analyses from multiple comparisons.
		- `roadmap_epigenomics_index.txt`: The index file of histone modifications needed by `chromatinPlots.R`.
	- `panTissue`: The scripts and annotation files for the Down Syndrome pan-tissue meta-analysis that utilizes [GAT](https://github.com/AndreasHeger/gat).
		- `annotations`: Folder with bed files for testing. References are in the manuscript.
		- `both`: Folder with GAT script and plotting script for each individual comparison.
		- `hyper_hypo`: Folder with GAT script and plotting script for each individual comparison split by hypermethylation and hypomethylation.
		- `panTissuePlot.R`: Plot of the analyses results for all comparisons from `runs_DS_GAT.sh`.
		- `run_DS_GAT.sh`: Main script for parallel runs of the Down syndrome pan-tissue meta-analysis for all comparisons. 

- `Figure 5`: RUNX1 profile.
	- `blockPlot`
		- `consensusBlocks.sh`: Create a consensus block profile from multiple comparisons for `RUNX1_blockPlot.R`.
		- `RUNX1_blockPlot.R`: Plot the RUNX1 consensus block using all samples.
	- `ChIP-seq_integration`:
		- `mmc3.xlsx`: RUNX1 ChIP-seq peaks. Reference is in the manuscript.
		- `RUNX1_ChIP_Overlap.R`: Analysis and plot of integration of WGBS data with RUNX1 ChIP-seq data.
	- `Transcription_Factor_Analysis`
		- `homer.sh`: Parallel HOMER known transcription factor motif analyses for multiple comparisons.
		- `marge.R`: Parsing and summary heatmap for HOMER analyses of multiple comparisons using MARGE.
		
## Reference

If you modify this code and utilize it please cite the publication:

Laufer BI, Hwang H, Jianu JM, Mordaunt CE, Korf IF, Hertz-Picciotto I, LaSalle JM. Low-Pass Whole Genome Bisulfite Sequencing of Neonatal Dried Blood Spots Identifies a Role for RUNX1 in Down Syndrome DNA Methylation Profiles. *Human Molecular Genetics*, 2020. **doi**: [10.1093/hmg/ddaa218](https://doi.org/10.1093/hmg/ddaa218)

## Acknowledgements

The development of these scripts was suppourted by a Canadian Institutes of Health Research (CIHR) postdoctoral fellowship [MFE-146824] and a [CIHR Banting postdoctoral fellowship](https://banting.fellowships-bourses.gc.ca/en/2018-2019-eng.html) [BPF-162684].  We would also like to thank [Matt Settles](https://github.com/msettles) for invaluable discussions related to some of the bioinformatic approaches utilized in this repository.
