This repository contains the code used for the analyses presented in the manuscript "BANDITS: Bayesian differential splicing accounting for sample-to-sample variability and mapping uncertainty", which presents the differential splicing R/Bioconductor package BANDITS.

Link to manuscript (currently on Biorxiv): https://doi.org/10.1101/750018

Link to BANDITS Bioconductor package: https://bioconductor.org/packages/release/bioc/html/BANDITS.html

The code is organized in four folders:

- `Simulated data 6 vs 6` contains the code to simulate the "6 vs 6" data used in the paper.
The simulated RNA-seq reads are available at FigShare (https://figshare.com/projects/RNA-seq_simulated_data_for_differential_transcript_usage_DTU_analyses/66275) with DOI `10.6084/m9.figshare.9467144`, `10.6084/m9.figshare.9692429` and `10.6084/m9.figshare.9692918`.

- `Alignment and quantification` contains the code to align and quantify reads.

- `Methods` contains the code to run all differential splicing methods.

- `Plots and Tables` contains the code to make all Figures and Tables available in the manuscript and its Supplementary. 
