# r-DMFA protocol which can be used to reproduce the results in the paper.

To reproduce the results presented in the paper, 
we have incorporated a supplementary module into the original computational framework. 
This module encompasses necessary modifications to enhance the analysis capabilities. 
We have also employed a customized version of the MATLAB implementation originally developed by 
Leighty and Antoniewicz in 2011, which is referenced in our paper. 
This dynamic model, along with the complete set of scripts needed for generating the study's results, 
has been compiled and made publicly accessible.

For researchers and practitioners interested in replicating or extending our study, 
the dynamic metabolic flux analysis (r-DMFA) protocol, including all the requisite scripts, 
can be downloaded from the following GitHub repository: https://github.com/ssbio/r-DMFA.

Please note that this repository contains all the updated scripts 
and the dynamic model necessary for the simulations discussed in the paper. 
By following the instructions and using the provided scripts, 
users can validate the findings or further explore the model's capabilities..


## Prerequisites

Before running the scripts, make sure to install all required dependencies listed below:
- IBM CPLEX solver - An academic license is obtainable for qualified institutions.
- Baron Solver interfaced with Matlab - MATLAB-BARON interface baron function available at [MINLP](https://www.minlp.com/mat-bar-baron). Licensed for a fee.
- Python 3
- os (standard library, no need to install)
- pandas
- sklearn (specifically sklearn.manifold for TSNE and sklearn.decomposition for PCA)
- matplotlib (specifically matplotlib.pyplot and mpl_toolkits.mplot3d.Axes3D)
- numpy
- seaborn
- adjustText (can be installed via conda from conda-forge or via pip)
- scipy (specifically scipy.cluster.hierarchy for linkage and dendrogram)
- sklearn.cluster for DBSCAN
- sklearn.metrics for silhouette_score
- MATLAB: Version 2016b or newer is required. 
The code was implemented using MATLAB 2022b on the University of Nebraska-Lincoln (UNL) Holland Computing Center’s (HCC) high-performance cluster, SWAN.

## Running the Workflow

This workflow consists of a series of scripts (./code) that should be run in sequential order from `S1` to `S7`. Below is the order of execution along with a brief description of each script:

### S1a_GetN15AbundanceProfiles.m
This script is responsible for obtaining N15 abundance profiles. It is the starting point of the analysis.

### S2a_Estimate_kcat_values.m
Estimates kcat values based on the proposed model in the paper and within the system.

### S3_GetWTFluxesAndConc.m
Calculates the wild type fluxes and concentrations, setting a baseline for further analysis.

### S4a_GetEffectiveEnzymeDemand.m
Determines the effective enzyme demand in the system, which can inform on system bottlenecks.

### S5_r_dmfa_WithFluxSampling.m
Performs flux sampling for robustness analysis of the dynamic metabolic flux analysis (DMFA) model.
work with small sample size, say 1E2 or as much as your computing resource can handle. 
We had more resources at our disposal and thus worked with very big file (4E6 dynamic flux profile samples...)

### S5b_GetResidualsCorrelationFromDynamicNFluxSamples.m
Analyzes residuals and their correlation derived from the DMFA, crucial for model validation.

### S5c_GetOutliersIndicesBasedonRhoValues.m
Identifies outlier indices based on the Rho statistic, aiding in quality control.

### S6a_GetFoldChangesOfDynFluxSamplesforClustering.m
Prepares fold changes of dynamic flux samples for clustering analysis.

### S6b_Cluster_FC_data.py
Clusters flux fold change data, providing insights into control patterns within the metabolic network for the samples size 1E6 as to reproduce Fig. 6 of the paper.

### S6c_Cluster_FC_data.py
A Python script that complements or iterates upon the clustering of flux fold change data for other sample size.

### S7a_GetMutantsFluxes_EnzKO.m
Analyzes the fluxes in various mutant strains with enzyme knockouts to study the effects on the metabolic network.

### S7b_GetMutantsFluxes_RxnsKO.m
Focuses on the flux changes due to reaction knockouts in mutant strains.

### S7c_GetZindexandRxnsCorrelationsforEnzsKO.m
Computes Z-index and correlates it with reactions, for identifying key regulatory reactions.

### S7d_GetZindexandRxnsCorrelationsforRxnsKO.m
An extension or variation of S7c script for deeper analysis at the reaction Level of Z-index and reaction correlations.

## Color Bar File

- `color_bar.id.mat`: This file is used for generating color bars in plots.

## Additional Scripts

- `DMFASampler_constr.m`: r-DMFA framework with constraints and sampler.
- `Get_EnrichmentProfiles.m`: Script for obtaining enrichment profiles.
- `knee_pt.m`: Determines the knee point which is used in various optimization and cutoff scenarios.
- `README`: This file, contains the necessary information to run the scripts.



The code has been developed and implemented with MATLAB version 2022b and Python version 3 installed. 
It is necessary to have near or suitable environments for the scripts to run properly. 
Additionally, the code is compatible only with MATLAB versions 2016b or newer.



## Support

For any questions or issues, please contact ssbio209@gmail.com


