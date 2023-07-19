# Gene Tree Estimation error in Lampyridae

This repository contains code to perform the analysis from the study: Höhna et al., Robustness of Divergence Time Estimation Despite Gene Tree Error: A Case Study of Fireflies (Coleoptera: Lampyridae).

## Data exploration

The first step of the analysis was to explore the data by computing the sequence coverage and the data summary statistics.

- rb scripts/0_extract_percentage_missing.Rev
- Rscript scripts/0_plot_percentage_missing.R
- rb scripts/0_extract_data_summary.Rev


## Analyses of single loci

We estimated the phylogeny for each locus separately. Therefore, we assumed a GTR+Gamma substitution model and an unrooted tree prior with exponential branch length distribution. We ran four MCMC replicates and performed a posterior predictive simulation. This all can be done automatically when you run the script *1_mcmc_single_gene_unrooted.Rev*. However, note that you need to prespecify several variables in RevBayes first, for example:
- GENE_NAME="Lampypridae_1"
- DATA_DIR="data_AHE"
- OUTPUT_DIR="1_unrooted_gene_trees/output"
- N_REPS=4
- NUM_MCMC_ITER=50000
- N_CHAINS=1

Afterwards, you should check for convergence using *Rscript scripts/check_convergence.R 4 1_unrooted_gene_trees/output 1_unrooted_gene_trees/results Lampypridae_1*.

Next, if the run has converged you should run the script `1_mcmc_summary.Rev` with the same arguments as before in RevBayes, which computes the posterior probabilities of the clades. Then, you should run `1_mcmc_summary.Rev` with the same arguments as before in RevBayes, which computes the summary statistics of the posterior predictive simulations.

Overall, you should repeat this for each locus if you found good settings, e.g., adequate length of MCMC runs.

## Divergence time estimation
In this set of analyses we estimate the divergence times using a concatenated analysis of a subset of the loci. To run this analysis, you need to specify first the following variables in RevBayes:
- DATA_DIR="data_AHE"
- OUTPUT_DIR="2_divergence_times/output"
- FOSSILS="photinus" (can be either 'photinus' or 'no_photinus')
- REP=1 (specify the replicate directly here, we used this for better parallelization)
- N_REPS=1 (since we specify the replicates directly, don't change this)
- NUM_MCMC_ITER=50000
- N_CHAINS=1
- DATA_FILTER="geneshopping" (can be either 'PP', 'low_GC_var', 'geneshopping' or 'high_coverage')
Then run the RevBayes script *2_mcmc_concatenated_rooted.Rev*. We recommend to use the parallel version of RevBayes for this. This analysis can take a long time and might have convergence issues. Be careful.


## Simulation study to assess gene tree estimation error
We performed a simulation study to see the gene tree estimation error on simulated data compared to empirical data. You can simulate gene tree in RevBayes using the script `3_sim_gene_trees.Rev` where you need to specify the variable `TRUE_POP_SIZE=1E5` or set it to any appropriate value. Note that we assume the species tree to be in millions of years, thus we internally also rescale the population size. Once you have the gene trees, you can simulate alignments for the gene tree using the script `3_sim_alignments.Rev`. Here, you need to specify before the following variables:
- TRUE_POP_SIZE <- 1E5
- NUM_GENES <- 436
- SIM_DIR <- "3_simulation_study_gene_trees/sim_trees"
- SIM_DIR <- "3_simulation_study_gene_trees/sim_alignments"
- SEED <- 12345
- RATE_VAR = "UCE"  (can be either 'UCE' or anything else)
Once you have the simulated alignments, you can run the script *1_mcmc_single_gene_unrooted.Rev* as above on each of the simulated datasets.
