#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RESOS  |  README FILE | December 2021
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PROJECT DESCRIPTION: 
This directory contains all the code from testing the RESOS. 

## FILES IN THE DIRECTORY

1. Functions
	1.1 DynamicPercentileSampling          -> Finds the optimal number of solutions  of a sampling based on a dynamic percentile computed with the difference 
		between the furthest and most centric solution per percentile. The inflection point of the curve is the optimal percentile. 
	1.2 EuclideanDistancesMaxMin           -> Computes the distance of one solution against every other solution of a sampling. It returns the furthest and central solution
		of the densest space of solutions computed with the dynamic percentile.  
	1.3 EuclideanDistancesSmallModels      -> Computes the distance of one solution against every other solution of a sampling. It returns the densest space of solutions 
		computed with the dynamic percentile. 
	1.4 FrequencyHistogramsMaxMin          -> Computes the frequency of each flux value in a solution and compares the total frequency value of a solution solution against every other solution of a sampling.
 			It returns the furthest and central solution of the densest space of solutions computed with the dynamic percentile.
	1.5 FrequencyHistogramsSmallModels     -> Computes the frequency of each flux value in a solution and compares the total frequency value of a solution solution against every other solution of a sampling.
 			It returns the densest space of solutions 
		computed with the dynamic percentile.


2. Models

	2.1 IMM_iCHO_R_33 -> Model used for the examples.  

3. Examples
	3.1 example_1     -> Applying RESOS algorithm with gpSampler
	3.2 example_2     -> Applying RESOS algorithm with ADBS 

