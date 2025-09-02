# pangenome-simulation

This script simulates gene presence/absence patterns across a phylogeny using a Wright–Fisher forward process with lineage-specific heterogeneity. A random tree is generated, and each gene is assigned a base fitness at the root. Genes are then seeded at the root, and their gain, loss, selection, and drift are simulated along branches. Branch-specific variation is introduced by an Ornstein–Uhlenbeck trait that modulates each gene’s effective selection coefficient.

The output includes a gene presence/absence matrix, a histogram of gene prevalence, and a phylogeny-ordered heatmap. A function is also provided to calculate pangenome fluidity. Example runs demonstrate “open” and “closed” pangenomes.
