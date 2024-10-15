This repository contains the codes for the paper: 
Zhao, Xingyuan, Ruyu Zhou and Fang Liu. "Enhancing Trade-offs in Privacy, Utility, and Computational Efficiency through MUltistage Sampling Technique (MUST)." arXiv preprint arXiv:2312.13389 (2023).

The codes for computing delta bounds given epsilon for the subsampled Gaussian mechanism over compositions are based on the 
repository https://github.com/DPBayes/PLD-Accountant, which contains the privacy profile analyses for Poisson subsampling and WOR with the Gaussian mechanism as the base mechanism but contains a bug in the case of WOR. We fixed the bug, allowing accurate privacy profile computations for WOR, and added the code for computing the privacy profile for WR and three MUST methods used with the Gaussian mechanism. 

1. The files bootstrap experiment.R, lm experiment.R, and logistic regression experiment.R contains the codes for utility experiments 1,2, and 3 respectively in Section 4.2. The code and results for the utility experiment of classification on the MNIST dataset are in MUST_MNIST.ipynb. 
2. The file privacy_amplification.R contains the R codes for the privacy amplification effect analyses and contour_3D_plot.R contains the codes for the contour and 3D plots. 
3. The FFT_revise_for_git.ipynb file contains the Python codes for privacy composition analyses and the adult_processed2.csv file contains the preprocessed Adult data used in utility experiment 3 logistic regression.



