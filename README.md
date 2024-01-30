Code used in generating the examples in the paper "Quantifying the impact of environmental changes on migratory species: a model perturbation framework".

MigModel.m is a MATLAB function file that computes the matrices that make up the annual cycle model for eight different migratory populations. 
  Input mm indicates the migratory model 
  Outputs are:  A (the seasonal matrix); Ahat (the seasonal survival matrix); n (the number of habitats/patches); c (the number of stage classes); s (the number of seasons); mP (contains all \mP_k the block matrices 									of the proportion of individuals migrating); D (D_{j,k} the demographic matrices); P (P^i_k the proportion of individuals migrtaing matrices); S (S^i_k the movement survival matrices); M (M^i_k = 										P^i_k.*S^i_k); mD (contains all \mD_k the block matrices of demography in all habitats); mS (contains all \mS_k the block matrices of movement survival); mM (contains all \mM_k).

PathwayMetrics.m is a MATLAB function file 
