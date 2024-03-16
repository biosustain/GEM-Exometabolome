#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ESSENTIAL REACTIONS & SYNTHETIC LETHAS - METHOD VALIDATION  |  README FILE | June 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DIRECTORY DESCRIPTION: 
This directory contains all the data/code & results from testing the RESOS method to identify essential reactions/synthetci lethals in the Ecoli GEM iJO1366. 

## FILES IN THE DIRECTORY
1. ExperimentalData (Directory)
	1.1 DATA/RAW (Directory) 
		Cotains the original list of essential genes and SL genes retrieved  from literature and the filtered list of only essential genes (LethalsGenesEco_Exp.csv). 
	1.2 DATA/RESULTS (Directory)
		Contains the essential reactions and SL reactions mapped into the GEM from the experimental essential genes/SL genes. 
	1.3 LethalsFromGenes.xml -> Script to get all posible ERs and SLs from the list of experimental genes.
	1.4 SynthicLethalsFromGenePairs.R -> Analysis of experimental SL.
	1.5 EssentialRxnsFromEssentialGenes_Exp.R -> Analysis of experimental essential reactions.

2. Simulations (Directory) 
	2.1 DATA/Results -> Contains the essential reactions/SL reactions obtained from each simulation. 
	2.2 DATA/Transcriptomic_Eco.csv -> Gene expression data used to constrain the model.
	2.3 DATA/iJO1366.mat.gz -> E coli GEM used in the simulations.
	2.4 AnalysisOfSimulations.R-> SynthicLethalsFromGenePairs.R -> Analysis of experimental SL.
	2.5 EcoliSim_Validations.m -> MATLAB script to perform simulations.
