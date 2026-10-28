# Coexistence rules for small, antagonistically interacting microbial communities

----

**Abstract.** 
The coexistence of diverse microbial communities despite the common presence of antimicrobial weapons presents a fundamental puzzle in ecology. 
To address this issue, we investigate the role of antibiotic-mediated interactions in driving microbial diversity using methods from graph theory and theoretical ecology. 
Our exhaustive analysis of small interaction graphs involving antibiotic production, resistance, and degradation reveals that two factors together determine whether an interaction pattern can support coexistence or not: a certain producer-sensitive-degrader (PSD) motif is critical, and a form of cyclicity between the action of different antibiotics is necessary.
Using individual-based simulations, we also explore the role of the PSD motif in spatially structured populations and demonstrate that community coexistence is robust over a wide range of antibiotic and degrader diffusivities. 
Our findings provide a deeper understanding of the interaction patterns that drive diversity in complex microbial communities.
Specifically, we emphasize how antagonism does not preclude biodiversity. 
These results offer clear pathways for cultivating synthetic microbial consortia, enabling the design of more effective strategies for manipulating microbial communities.

----

This repository contains the online material, codes for generating figures and the relevant simulations for reproducing the work with the associated manuscript.

The structure of the repository is as follows: 

./FigureSources - final Omnigraffle files (and some auxiliary files) for making each of the figures in the main and supplementary texts.

./SimulationCode - scripts for analysis and creation of auxiliary files for figures. This has two subfolders that are explained within.

Each subfolder contains its own README file to help the reader navigate.

----
## Additional analytical note

[Cyclic-pair necessity in the three-strain full-PSD class](notes/full-psd-cyclicity/README.md) gives a proposed proof and reproducible checks for arbitrary finite channel multiplicities and positive effective growth. It concerns locally asymptotically stable interior fixed points of the well-mixed map; the broader cyclicity conjecture remains unresolved. External author review is pending.
