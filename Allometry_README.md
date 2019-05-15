**Paper authors: A Shankar, DR Powers, LM Davalos, CH Graham**

Code by: A Shankar, github/nushiamme; contact: nushiamme\<at\>gmail\<dot\>com for questions about code/datasets

#### Code organisation

-   **Allometry_phylo.R** - Needs input files *"DLW_TableS1.csv"* and *"hum294.tre"* (contact A Shankar for access); contains code plotting and analysing allometric relationship in hummingbirds, accounting for phylogenetic relatedness.
    -   *Table 1*: Results from PGLS models (Brownian motion and Ornstein-Uhlenbeck models) of log(Daily energy expenditure) vs log(mass)
    -   *Table 2*: Results of the MCMCglmm models with the Brownian-motion tree and without a phylogenetic tree. All the models used individual values, except one model which used species means.
    -   *Fig. 1*: A log-log plot of individual and species-level hummingbird daily energy expenditure (kJ) vs. mass (g), including values from this study as well as from the literature. Includes species means and individual points.

#### Packages you will need:

    + MCMCglmm
    + nlme
    + ape
    + geiger
    + caper
    + phytools
    + RColorBrewer
    + ggplot2