# pangenome-simulation

This script simulates gene presence/absence patterns across a phylogeny using a Wright–Fisher forward process with lineage-specific heterogeneity. A random tree is generated, genes are seeded at the root, and their gain, loss, selection, and drift are simulated along branches. Branch-specific variation is introduced by an Ornstein–Uhlenbeck trait that modulates gain and loss rates.

The output includes a gene presence/absence matrix, a histogram of gene prevalence, and a phylogeny-ordered heatmap. A function is also provided to calculate pangenome fluidity. Example runs demonstrate “open” and “closed” pangenomes.

## Future to do

Explore different ways of modelling lineage heterogeneity. For example, give each gene a base fitness at the root, then each lineage modulates that fitness randomly.
