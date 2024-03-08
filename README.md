Code used in generating the examples in the paper "Quantifying the impact of environmental changes on migratory species: a model perturbation framework".

MigModel.m is a MATLAB function file that computes the matrices that make up the annual cycle model for eight different migratory populations. 
  Input mm indicates the migratory model 
  Outputs are:  A (the seasonal matrix); Ahat (the seasonal survival matrix); n (the number of habitats/patches); c (the number of stage classes); s (the number of seasons); mP (contains all \mP_k the block matrices of the proportion of individuals migrating); D (D_{j,k} the demographic matrices); P (P^i_k the proportion of individuals migrtaing matrices); S (S^i_k the movement survival matrices); M (M^i_k = P^i_k.*S^i_k); mD (contains all \mD_k the block matrices of demography in all habitats); mS (contains all \mS_k the block matrices of movement survival); mM (contains all \mM_k).

PathwayMetrics.m is a MATLAB function file that computes the subpopulation pathway contribution metrics (Cp) and the metapopulation pathway contribution metrics (Cptilde) for the pathways specified in PATH, see https://doi.org/10.1016/j.ecolmodel.2022.110056 for details about how they are computed. 
  Inputs are: c (the number of stage classes); n (the number of habitats/patches); s (the number of seasons); A (the seasonal matrix); Ahat (the seasonal survival matrix); Phi (stores the seasons in which a path is specified); PATH (stores the pathways of interest); mP (contains all \mP_k the block matrices of the proportion of individuals migrating).
  Outputs are: Cp (the subpopulation pathway contribution metrics); Cptilde (the metapopulation pathway contribution metrics); AAhat (the product of the A and Ahat matrices for the specified path); AAhatP (the product of the A, Ahat and P matrices for the specified path).

DistinctPaths.m is a MATLAB function file that computes a matrix containing all the distinct paths for a model with s seasons in the annual cycle and n habitats.
  Inputs are: n (the number of habitats/patches); s (the number of seasons); sp (the seasons that have paths specified).
  Output PATH is a matrix whose rows contain the distinct pathways.

These three MATLAB function files are used in the MATLAB scripts that generate the examples in the paper. 

Example 3.1: Simple Hypothetical Model Perturbation Framework: 
Use file pertCmetrics.m to compute the change in contribution metrics using the perturbation framework. 

Example 3.2: Simple Hypothetical Model Sensitivity Framework:
Use file sensCmetrics.m to compute sensitivities and elasticities of the contribution metrics.

Example 3.3: Monarch Butterfly Model Sensitivity of Life Rates in Season 6:
Use file MonarchExamplePPCM.m to compute sensitivities and elasticities of the contribution metrics.
Use file sensLambdaMonarch.m to compute sensitivities and elasticities of the asymptotic growth rate. 

Example 3.4: Monarch Butterfly Model Threats and Proposed Management Actions:
Use file MonarchExamplePPCM2.m to compute change in contribution metrics following threats and threats and conservation. Contains code to plot figure using matlab. Also contains code to organise the data for plotting in R. 
Use file Monarch2Plot.R to plot figure as seen in paper. 
