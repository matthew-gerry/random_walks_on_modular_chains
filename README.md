This repository contains a set of Matlab functions and scripts that can be used to produce all the figures found in the article 'Random Walks on Modular Chains: Detecting Structure through Statistics' by Matthew Gerry and Dvira Segal, Phys. Rev. E 108, 024135 (2023) (available: https://arxiv.org/abs/2305.09759 and https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.024135). This includes the implementation of full counting statistics for a random walk on a chain with differing site-site transition rates that repeat in a modular fashion, amounting to a spatial periodicity in the structure. In addition, some of the Matlab files included here numerically generate probability distributions over the sites of such random walks to aid in providing inuition for why the particular analytic results obtain. The goal of this project has been to understand how studying the statistics of a Markov jump process can play an important role in determining its underlying structure.

All of the figures are generated by running the scripts cumulant_plots_modRW.m and pdf_plots_modRW.m, and users may edit the Matlab scripts to change the values of physical parameters or even the structural characteristics of the chain on which the random walk, which is the subject of this study, plays out.