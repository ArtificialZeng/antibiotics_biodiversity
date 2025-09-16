This folder contains the scripts necessary for analysis of the mixed-inhibition zone model. 

SCRIPTS:

`get_NRcommunities.m` can be used to get all non-redundant communities with N species and M antibiotics, which can be defined within the script. This will contain communities with repeated rows and columns, these are filtered out in `SearchSpaceOfCommunities.m` below. This is not the most efficient way of doing this, but it works.

`SearchSpaceOfCommunities.m` explores the space of communities found above and saves a .mat file containing the indices (one time for each fixed point found by the solver) in Results/Search-Space/. The communities that these indices correspond to can be found by calling these indices in the appropriate (N,M) file from ./NRcommunities.

`heatmap_basin_size.m`: for a given input community and parameters, calculate the size of the basin of attraction for its fixed point (s, if it has more than one) for different cp/cr and cp/cd values.

`heatmap_spectralradius.m`: for a given input community and parameters, calculate the spectral radius of the Jacobian evaluated at a fixed point (if it exists) for different cp/cr and cp/cd values.

`heatmap_spectralradius_KPKD.m`: for a given input community and parameters, calculate the spectral radius of the Jacobian evaluated at a fixed point (if it exists) for different KP and KD values.

`heatmap_final_PSR.m` and `heatmap_final_PSD.m`: calculate the final state of the trajectories for the communities [P,S,R] and [P,S,D] for different cost parameters. 

FOLDERS:

./functions contains the functions (written as individual .m files) necessary to run the above scripts.

./NRcommunities contains files that store non-redundant lists of communities for a few N,M value pairs.

./Results contains heatmaps etc. that we use to make the figures.


