#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ESSENTIAL REACTIONS & SYNTHETIC LETHAS - METHOD VALIDATION  |  README FILE | June 2022
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PROJECT DESCRIPTION: 
This directory contains all the data/code & results from testing the ROOM+Sampling method to identify essential reactions/synthetci lethals. 

## FILES IN THE DIRECTORY

1. ExperimentalData (Directory)

	1.1 RAW (Directory) 
		1.1.1 KEIO_genes.csv -> list of all the genes in the KEIO collection. - This must be filter to only the essential genes.
 		1.1.2 SLGenes.csv -> list of all experiments performed to test SLs - This must be filtered to the top 5 percentile of the log(Q/R)
	1.2 SingleLethalsGenes_Exp_7Jun.csv -> List of essential genes. 
	1.3 SLGenes_Exp_7Jun.csv -> List of SL genes.
	1.4 LethalsFromGenes.xml -> Gets all posible ERs and SLs from the list of genes (1.2 and 1.3)
	1.5 tableFinalGene1.csv -> List of all rxns associated with gene 1 of the SL genes - This must be filtered using GPRs and combined with tableFinalGene2. 
	1.6 tableFinalGene2.csv -> List of all rxns associated with gene 2 of the SL genes - This must be filtered using GPRs and combined with tableFinalGene1.
 	1.7 PotentialLethalRxns_7Jun.csv -> List of all rxns associated with gene essential genes - This must be filtered using GPRs. 
	1.8 EssentialRxnsFromEssentialGenes_Exp.R ->Gets the experimental essential reactions using GPRs - need the matlab list of rxns (1.7) and the list of genes (1.2) 	
	1.9 SynthicLethalsFromGenePairs -> Gets the list of SLs using GPRs - need the matlab lists (1.5,1.6) and the lists of genes (1.2,1.3)
	1.10 LethalRxns_Exp_7Jun -> FINAL list of lethal reactions - use this to compare the results. 
	1.11 Synthetic lethalRxns_Exp_7Jun -> FINAL list of SLs - use this to compare results. 


2. ROOM+Compact+Loopless (Directory) 

	2.1 EcoLooplessCompacted_28Jun.mat -> Matlab workspace for compacted+loopless+ROOM test - the models and final results are here.
	2.2 Analysis_28Jun.R -> Comparison analysis between exp data and Compacted+loopless+ROOM and between exp data and OriginalModel+FastSL
	2.3 SL1_CompactedLoopless_28Jun.csv -> List of rxns 1 of 2 for SLs from Compacted+Loopless+ROOM - Must be combined with SL2_SL2_CompactedLoopless_28Jun.csv 
	2.4 SL2_CompactedLoopless_28Jun.csv -> List of rxns 2 of 2 for SLs from Compacted+Loopless+ROOM - Must be combined with SL1_CompactedLoopless_28Jun.csv 
	2.5 ERs_Compacted_Loopless_28Jun.csv -> List of essential reactions from Compacted+Loopless+ROOM test. 
	2.6 SLs_Original_28Jun.csv -> List of SLs from originalModel+FastSL test. 
	2.7 ERs_Original_28Jun.csv -> List of ERs from the OriginalModel+FastSL test. 
	2.8 FinalPipelineEcoli_ROOM.m -> Loopless+Compacted+ROOM test (code) - 2.1 is the workspace at the end of this. - The results of ONLY the parallel SL are results  of the OriginalModel+FastSL. 

3. FastSL+Loopless/Loopless+Compact (Directory)

	3.1 Compacted_or_Loopless+FastSL_analysis.R -> Comparison analysis between exp data vs Compacted+Loopless+FastSL and Loopless+FastSL vs exp data.
	
	3.2 SLs_Loopless_28Jun.csv -> SL from Loopless+FastSL test. 
	3.3 ERs_Loopless_28Jun.csv -> ER from Loopless+FastSL test. 
	3.4 SL1_Compacted_28Jun -> Reactions 1 of 2 from SL from Compacted+Loopless+FastSL test - must be combines with 3.5 
	3.5 SL2_Compacted_28Jun -> Reactions 2 of 2 from SL from Compacted+Loopless+FastSL test - must be combined with 3.4 
	3.6 ERs_Compacted_28Jun -> ERs from Compacted+Loopless+FastSL test. 	
 	3.7 modelLoopless.mat -> Loopless model used for the Loopless+FastSl test
	3.8 model_compacted.mat -> Loopless+Compacted model used for the Loopless+Compacted+FastSl test. 
	3.9 OtherTest.m -> script for Loopless+Compacted+FastSl test and Loopless+FastSl test.

4.Test5 (Directory)
	4.1 ERs_ROOMALL_28Jun.csv -> ERs from Loopless+Compacted+ROOM(all rxns)
	4.2 Test5_Results_28Jun.csv -> Results analysis 
	4.3 Test5.m -> MATLAB script for Loopless+Compacted+ROOM(all rxns)

5. metabolicEco_Average.csv -> Trancriptomic data used for CiMAT/iMAT. 
6. reportPlots.r -> Plots of the results (Not really important) 
7. EcoTest-29Jun.xlsx -> Contingency tables for all the tests. 


