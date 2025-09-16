Pipeline for the creation of each of the figures. Figures are made with the software OmniGraffle.

There is a .graffle file for each figure, which is the main method of compiling sub-components of the figure and making aesthetic changes.

--- 

Main text figures:

Figure 1: only conceptual. Everything required for this figure is in the graffle and associated Mathematica notebook of the same name.

Figure 2: Panels (a.1-a.4) and (b.1-b.4) are made by running `SimulationCode/mixed_inhibition_zone/heatmap_basin_size.m` for appropriate community and parameters. The output is saved to the same folder. Results are stored in `SimulationCode/mixed_inhibition_zone/Results/2antibiotic_attractorpresence`.

Figure 3: Only conceptual. Results are obtained by running `SimulationCode/mixed_inhibition_zone/SearchSpaceOfCommunities.m` with N=M=3 and appropriate parameters. Results are stored in `SimulationCode/mixed_inhibition_zone/Results/Search-Space`.

Figure 4: Panels (a.1, b.1, c.1) are made by running `SimulationCode/mixed_inhibition_zone/heatmap_spectralradius.m` for appropriate community and parameters. Panels (a.2, b.2, c.2) are made by running `SimulationCode/mixed_inhibition_zone/heatmap_basin_size.m` for appropriate community and parameters. Results are stored in `SimulationCode/mixed_inhibition_zone/Results/cyclicity-complexity`.

Figure 5: Simulation model is in `SimulationCode/spatial_model/HPC_spatial_model.py`, and data is analysed (plots are made) using the script `SimulationCode/spatial_model/analyse_spatial_model.py`.   

Figure 6: the heatmap on the left is again made using `analyse_spatial_model.py`,and the simplex on the right is colour_simplex.py which is made using the script colour_simplex.py

--- Supplementary figures

Figures A1-A4 are only conceptual and standalone. Everything required to make the figure can be found in the corresponding .graffle file.

Figure A5: Top row of heatmaps made using `SimulationCode/mixed_inhibition_zone/heatmap_final_PSR.m` for appropriate parameters; bottom row using `SimulationCode/mixed_inhibition_zone/heatmap_final_PSD.m`. Results are stored in `SimulationCode/mixed_inhibition_zone/Results/1antibiotic-PSD_PSR`.

Figure A6: Panel (b) made by running `SimulationCode/mixed_inhibition_zone/heatmap_spectral.m` for appropriate community and parameters. Panel (a) using `SimulationCode/mixed_inhibition_zone/heatmap_spectralradius_KPKD.m`. Results are stored in `SimulationCode/mixed_inhibition_zone/Results/cyclicity-complexity`.

Figure A7: Trajectories are obtained by running `SimulationCode/mixed_inhibition_zone/specific_communities.m` with appropriate community and parameters. Results are stored in `SimulationCode/mixed_inhibition_zone/Results/trajectories`.


