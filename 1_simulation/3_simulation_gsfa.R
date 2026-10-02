# GSFA simulation: generate data with PerturbVI, fit GSFA, report metrics.

# --- scenarios to run ---
#   b:      0.05, 0.10, 0.15, 0.20
#   g:      50, 100, 150, 200
#   z:      3, 4, 5, 6
#   l:      200, 400, 600, 800
#   mis_g:  0.1, 0.2, 0.3, 0.4, 0.5
scenarios <- list(
  c("b", 0.05), c("b", 0.10), c("b", 0.15), c("b", 0.20),
  c("g", 50),   c("g", 100),  c("g", 150),  c("g", 200),
  c("z", 3),    c("z", 4),    c("z", 5),    c("z", 6),
  c("l", 200),  c("l", 400),  c("l", 600),  c("l", 800),
  c("mis_g", 0.1), c("mis_g", 0.2), c("mis_g", 0.3), c("mis_g", 0.4), c("mis_g", 0.5)
)
n_sim <- 80
seed_start <- 0       
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
out_dir <- file.path(dirname(if (length(script)) script else "."), "results")

library(reticulate)
library(GSFA)

use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)

py_run_string('
import numpy as np
from jax import random
import perturbvi
from procrustes import orthogonal

def procrustes_error(A, B):
    return orthogonal(np.asarray(A), np.asarray(B), scale=True, pad=True).error

def compute_sen_spec(overall, lfsr):
    overall = np.asarray(overall)
    lfsr = np.asarray(lfsr)
    nz = int((overall != 0).sum())
    z = int((overall == 0).sum())
    sel_nz = int(np.logical_and(overall != 0, lfsr < 0.05).sum())
    sel_z = int(np.logical_and(overall == 0, lfsr < 0.05).sum())
    return sel_nz / nz, sel_z / z
')

run_one <- function(sweep, value) {
  # data-generating constants
  n_dim <- 3000L
  p_dim <- 4000L
  data_z_dim <- 4L
  data_g_dim <- 100L
  data_l_dim <- 150L
  data_b <- 0.20
  fit_k <- 4L

  if (sweep == "b") {
    data_b <- value
    out_name <- paste0("gsfa_b_", value, ".csv")
  } else if (sweep == "g") {
    data_g_dim <- as.integer(value)
    out_name <- paste0("gsfa_g_", as.integer(value), ".csv")
  } else if (sweep == "z") {
    fit_k <- as.integer(value)
    out_name <- paste0("gsfa_z_", as.integer(value), ".csv")
  } else if (sweep == "l") {
    data_l_dim <- as.integer(value)
    out_name <- paste0("gsfa_l_", as.integer(value), ".csv")
  } else if (sweep == "mis_g") {
    remove_prop <- value
    out_name <- paste0("gsfa_misg_", as.integer(value * 100), ".csv")
  } else {
    stop("unknown sweep: ", sweep)
  }

  if (sweep == "mis_g") {
    niter <- 200L
    used_niter <- 100L
  } else {
    niter <- 2000L
    used_niter <- 1000L
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  gen_code <- function(seed) {
    sprintf("
Z, W, X, G, beta = perturbvi.generate_sim(
    seed=%d, l_dim=%d, n_dim=%d, p_dim=%d, z_dim=%d, g_dim=%d, b_sparsity=%.6f)
X = np.asarray(X - X.mean(axis=0))
G = np.asarray(G)
W = np.asarray(W)
Z = np.asarray(Z)
beta = np.asarray(beta)
", seed, data_l_dim, n_dim, p_dim, data_z_dim, data_g_dim, data_b)
  }

  results <- list()
  for (i in seq_len(n_sim)) {
    seed <- seed_start + (i - 1L)
    py_run_string(gen_code(seed))

    if (sweep == "mis_g") {
      py_run_string(sprintf("
remove_idx = np.asarray(random.choice(random.PRNGKey(%d), %d, shape=(%d,), replace=False))
G_sub = np.delete(G, remove_idx, axis=1)
beta_sub = np.delete(beta, remove_idx, axis=0)
overall_sub = (beta_sub @ W).T
", seed, data_g_dim, as.integer(data_g_dim * remove_prop)))

      fit_correct <- fit_gsfa_multivar(Y = py$X, G = py$G, K = fit_k,
        prior_type = "mixture_normal", init.method = "svd",
        prior_w_s = 50, prior_w_r = 0.2, prior_beta_s = 5, prior_beta_r = 0.2,
        niter = niter, used_niter = used_niter, verbose = FALSE, return_samples = FALSE)
      fit_miss <- fit_gsfa_multivar(Y = py$X, G = py$G_sub, K = fit_k,
        prior_type = "mixture_normal", init.method = "svd",
        prior_w_s = 50, prior_w_r = 0.2, prior_beta_s = 5, prior_beta_r = 0.2,
        niter = niter, used_niter = used_niter, verbose = FALSE, return_samples = FALSE)

      lfsr_c <- fit_correct$lfsr[, -ncol(fit_correct$lfsr), drop = FALSE]
      lfsr_m <- fit_miss$lfsr[, -ncol(fit_miss$lfsr), drop = FALSE]
      py$lfsr_c <- lfsr_c
      py$lfsr_m <- lfsr_m
      py_run_string("
lfsr_c_sub = np.delete(np.asarray(lfsr_c), remove_idx, axis=1)
sc = compute_sen_spec(overall_sub, lfsr_c_sub)
sm = compute_sen_spec(overall_sub, np.asarray(lfsr_m))
")
      results[[i]] <- list(sensitivity_correct = py$sc[[1]],
                           specificity_correct = py$sc[[2]],
                           sensitivity = py$sm[[1]],
                           specificity = py$sm[[2]])
    } else {
      fit <- fit_gsfa_multivar(Y = py$X, G = py$G, K = fit_k,
        prior_type = "mixture_normal", init.method = "svd",
        prior_w_s = 50, prior_w_r = 0.2, prior_beta_s = 5, prior_beta_r = 0.2,
        niter = niter, used_niter = used_niter, verbose = FALSE, return_samples = FALSE)

      pm <- fit$posterior_means
      beta_hat <- pm$beta_pm
      beta_hat <- beta_hat[1:(nrow(beta_hat) - 1), , drop = FALSE]
      z_hat <- pm$Z_pm
      w_hat <- pm$W_pm
      lfsr <- fit$lfsr[, -ncol(fit$lfsr), drop = FALSE]
      py_run_string("overall = (beta @ W).T")

      beta_err <- py$procrustes_error(beta_hat, py$beta)
      z_err <- py$procrustes_error(z_hat, py$Z)
      w_err <- py$procrustes_error(w_hat, t(py$W))
      roc <- py$compute_sen_spec(py$overall, lfsr)
      results[[i]] <- list(beta_err = beta_err, z_err = z_err, w_err = w_err,
                           sensitivity = roc[[1]], specificity = roc[[2]])
    }
    cat(sweep, value, "seed", seed, "done\n")
  }

  df <- do.call(rbind, lapply(results, function(x) as.data.frame(x)))
  write.csv(df, file.path(out_dir, out_name), row.names = FALSE)
  cat("wrote", file.path(out_dir, out_name), "\n")
}

for (s in scenarios) {
  run_one(s[1], as.numeric(s[2]))
}
