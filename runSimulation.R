# -----------------------------------------------------------------------------
# Simulating presence/absence of independent genes down a fixed phylogeny.
# -----------------------------------------------------------------------------
# - Propagation with Wright–Fisher forward simulation on every branch.
#
# - Per-generation selection (haploid, carriers fitness = 1 + s), gain/loss,
# then binomial sampling with Ne (drift). Parent state is treated as 0/1.
#
# - Rates of base_gain and base_loss are on the log scale; q = exp(base_*) gives
# instantaneous per-unit-time rates which are converted to per-generation probs.
#
# - Branch heterogeneity is simulated by determining a per-gene base fitness at 
# the root (s_root), then each lineage modulates that fitness via an 
# Ornstein-Uhlenbeck trait (sig2, alpha=1, theta=0).
#
# - Root seeding where each gene has root_prob chance of being present at the root.
# -----------------------------------------------------------------------------

library(ape)
library(phytools)
library(tidyverse)

set.seed(123)

# -------- WF branch propagation --------
propagate_branch <- function(p0, bl, q01_eff, q10_eff, s, Ne, gens_per_unit) {
  p <- min(max(as.numeric(p0), 0), 1)
  G <- max(1, round(bl * gens_per_unit))
  for (g in seq_len(G)) {
    if (p == 0) {
      p_sel <- 0
    } else if (p == 1) {
      p_sel <- 1
    } else {
      w1 <- 1 + s; w0 <- 1
      p_sel <- (p * w1) / (p * w1 + (1 - p) * w0)
    }
    p_mut <- p_sel * (1 - q10_eff) + (1 - p_sel) * q01_eff
    p_mut <- min(max(p_mut, 0), 1)
    if (is.infinite(Ne) || Ne <= 0) {
      p <- p_mut
    } else {
      k <- rbinom(1, Ne, p_mut); p <- k / Ne
    }
    if (p == 0 || p == 1) break
  }
  p
}

# -------- Simulate one gene --------
simulate_gene <- function(tree, params) {
  q01_inst <- exp(params$base_gain)
  q10_inst <- exp(params$base_loss)
  ntips_local <- length(tree$tip.label)
  states <- setNames(rep(NA, ntips_local), tree$tip.label)
  
  ou_s <- fastBM(tree, a = 0, sig2 = params$s_sig2, alpha = 1, theta = 0, internal = TRUE)
  
  root_p <- if (!is.null(params$root_prob)) params$root_prob else 0.5
  state_root <- rbinom(1, 1, root_p)
  
  traverse <- function(node, parent_state) {
    children <- tree$edge[tree$edge[,1] == node, 2]
    if (length(children) == 0) {
      states[tree$tip.label[node]] <<- parent_state
    } else {
      for (child in children) {
        bl <- tree$edge.length[which(tree$edge[,2] == child)]
        node_name <- if (child <= ntips_local) tree$tip.label[child] else as.character(child)
        
        trait_s <- ou_s[node_name]
        s_eff <- params$s_root + trait_s
        s_eff <- pmin(pmax(s_eff, -0.99), 10)
        
        q01_eff <- 1 - exp(-q01_inst / params$gens_per_unit)
        q10_eff <- 1 - exp(-q10_inst / params$gens_per_unit)
        q01_eff <- pmin(pmax(q01_eff, 0), 1); q10_eff <- pmin(pmax(q10_eff, 0), 1)
        
        p_parent <- if (parent_state == 1) 1 else 0
        
        p_child <- propagate_branch(
          p_parent, bl, q01_eff, q10_eff, s_eff, params$Ne, params$gens_per_unit
        )
        child_state <- rbinom(1, 1, p_child)
        traverse(child, child_state)
      }
    }
  }
  
  traverse(ntips_local + 1, state_root)
  list(states = states)
}

# -------- Run simulation --------
run_simulation <- function(params) {
  pa_mat <- matrix(0, nrow = length(params$tree$tip.label), ncol = params$G,
                   dimnames = list(params$tree$tip.label, paste0("gene", 1:params$G)))
  
  for (g in seq_len(params$G)) {
    params_gene <- params
    if (!is.null(params$s)) {
      params_gene$s_root <- params$s
    } else {
      s_mean <- if (!is.null(params$s_root_mean)) params$s_root_mean else 0
      s_sd   <- if (!is.null(params$s_root_sd)) params$s_root_sd else 0.01
      params_gene$s_root <- rnorm(1, mean = s_mean, sd = s_sd)
    }
    
    res <- simulate_gene(params$tree, params_gene)
    pa_mat[, g] <- res$states
  }
  
  prevalence <- colSums(pa_mat) / length(params$tree$tip.label)
  
  # ---- Histogram PDF ----
  prevalence_nonzero <- prevalence[prevalence > 0]
  hist_file <- paste0(params$tag, "_hist.pdf")
  pdf(hist_file, width = 7, height = 5)
  print(
    ggplot(tibble(prevalence = prevalence_nonzero), aes(x = prevalence)) +
      geom_histogram(binwidth = 0.02, fill = "gray70", color = "black") +
      labs(title = paste(params$tag, "prevalence"),
           x = "Fraction of tips carrying gene", y = "Count") +
      theme_minimal(base_size = 14)
  )
  dev.off()
  
  # ---- Phylogeny + heatmap PDF ----
  tree_file <- paste0(params$tag, "_tree.pdf")
  pdf(tree_file, width = 10, height = 8)
  gene_order <- order(colSums(pa_mat), decreasing = TRUE)
  phytools::phylo.heatmap(
    params$tree,
    pa_mat[, gene_order],
    standardize = FALSE,
    fsize = 0.6,
    colors = c("white", "black"),
    legend = FALSE,
    outline = FALSE
  )
  title(main = paste(params$tag, "tree and heatmap"), line = -1)
  dev.off()
  
  invisible(list(pa = pa_mat, prevalence = prevalence_nonzero,
                 hist_file = hist_file, tree_file = tree_file))
}

# -------- Pangenome fluidity --------

# Thanks to Anna Dewar for the function :-)

calc_pangenome_fluidity <- function(pa_mat) {
  fluidity_pair_eq <- function(x, y) {
    x <- as.integer(x); y <- as.integer(y)
    Mk <- sum(x); Ml <- sum(y)
    Uk <- sum(x == 1 & y == 0)
    Ul <- sum(x == 0 & y == 1)
    (Uk + Ul) / (Mk + Ml)
  }
  
  pair_combos <- combn(rownames(pa_mat), 2, simplify = FALSE)
  
  pair_vals <- sapply(pair_combos, function(p) {
    fluidity_pair_eq(pa_mat[p[1], ], pa_mat[p[2], ])
  })
  
  (2 / (nrow(pa_mat) * (nrow(pa_mat) - 1))) * sum(pair_vals)
}

# -------- Params explanation --------
# G - number of genes to simulate
# root_prob - probability gene is present at root
# gens_per_unit - generations per unit branch length
# Ne - effective population size
# base_gain - log-rate of gene gain per branch (instantaneous)
# base_loss - log-rate of gene loss per branch (instantaneous)
# s_root_mean / s_root_sd - per-gene base fitness (sampled per gene)
# s_sig2 - OU variance for the trait that modulates s along the tree
# (legacy) s - keep the same selection for all genes/branches

# -------- Examples --------
setwd("~/Desktop/pangenome-simulation")

ntips <- 100
tree <- rcoal(ntips)

# -------- Open pangenome example --------
params_open <- list(
  tag = "Open",
  tree = tree,
  G = 500,
  Ne = 1e9,
  gens_per_unit = 20,
  base_gain = -2,
  base_loss = -2,
  s_root_mean = 0,
  s_root_sd   = 0.005,
  s_sig2 = 0.01,
  root_prob = 0.5
)

sim_open <- run_simulation(params_open)
open_fluidity <- calc_pangenome_fluidity(sim_open$pa)

# -------- Closed pangenome example --------
params_closed <- list(
  tag = "Closed",
  tree = tree,
  G = 500,
  Ne = 1e9,
  gens_per_unit = 20,
  base_gain = -4,
  base_loss = -4,
  s_root_mean = 0,
  s_root_sd   = 0.005,
  s_sig2 = 0.01,
  root_prob = 0.5
)

sim_closed <- run_simulation(params_closed)
closed_fluidity <- calc_pangenome_fluidity(sim_closed$pa)
