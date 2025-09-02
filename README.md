# pangenome-simulation

This script simulates gene presence/absence patterns across a phylogeny using a Wright–Fisher forward process with lineage-specific heterogeneity. A random tree is generated, and each gene is assigned a base fitness at the root. Genes are then seeded at the root, and their gain, loss, selection, and drift are simulated along branches. Branch-specific variation is introduced by an Ornstein–Uhlenbeck trait that modulates each gene’s effective selection coefficient.

The output includes a gene presence/absence matrix, a histogram of gene prevalence, and a phylogeny-ordered heatmap. A function is also provided to calculate pangenome fluidity. Example runs demonstrate “open” and “closed” pangenomes.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fgithub.com%2Fwtmatlock%2Fpangenome-simulation%2F/main?urlpath=%2Fdoc%2Ftree%2FrunSimulation.R)

