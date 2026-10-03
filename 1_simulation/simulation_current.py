"""Simulate data and fit PerturbVI across the manuscript sweeps.

Scenarios (80 seeds, seeds 0..79):
  b_sparsity: 0.05, 0.10, 0.15, 0.20   (g=100, K=4, L=150)
  g_dim:      50, 100, 150, 200        (b=0.20, K=4, L=150)
  z_dim (K):  3, 4, 5, 6               (true K=4, b=0.20, g=100)
  l_dim (L):  100, 125, 150, 175, 200  (true L=150, b=0.20, g=100)

Usage:
  python simulation_current.py --b_sparsity 0.05 --out perturbvi_b_0.05.csv
  python simulation_current.py --g_dim 100 --out perturbvi_g_100.csv
  python simulation_current.py --z_dim 5 --out perturbvi_z_5.csv
  python simulation_current.py --l_dim 100 --out perturbvi_l_100.csv
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from jax import config
from jax.experimental import sparse
from procrustes import orthogonal

import perturbvi

config.update("jax_enable_x64", True)
config.update("jax_default_matmul_precision", "highest")

OUT_DIR = Path(__file__).resolve().parent / "results"


def compute_sen_spec(overall, lfsr):
    """Sensitivity/specificity of LFSR < 0.05 against nonzero entries of overall."""
    non_zero = (overall != 0).sum()
    zero = overall.size - non_zero
    selected_non = np.logical_and(overall != 0, lfsr < 0.05)
    selected_zero = np.logical_and(overall == 0, lfsr < 0.05)
    return selected_non.sum() / non_zero, selected_zero.sum() / zero


def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate and fit PerturbVI on one scenario."
    )
    p.add_argument("--n_dim", type=int, default=3000)
    p.add_argument("--p_dim", type=int, default=4000)
    p.add_argument("--true_l_dim", type=int, default=150,
                   help="single effects per factor used by the data generator")
    p.add_argument("--l_dim", type=int, default=150,
                   help="single effects per factor used when fitting")
    p.add_argument("--true_z_dim", type=int, default=4,
                   help="number of latent factors used by the data generator")
    p.add_argument("--z_dim", type=int, default=4,
                   help="number of latent factors used when fitting")
    p.add_argument("--g_dim", type=int, default=100)
    p.add_argument("--b_sparsity", type=float, default=0.20)
    p.add_argument("--p_prior", type=float, default=0.15)
    p.add_argument("--tau", type=float, default=1.0)
    p.add_argument("--tol", type=float, default=1e-3)
    p.add_argument("--max_iter", type=int, default=600)
    p.add_argument("--init", choices=["auto", "pca", "random"], default="auto")
    p.add_argument("--n_sims", type=int, default=80, help="number of seeds")
    p.add_argument("--seed_start", type=int, default=0)
    p.add_argument("--out", type=str, default="simulation_current.csv")
    return p.parse_args()


def main():
    args = parse_args()

    init = args.init
    if init == "auto":
        init = "pca" if args.z_dim <= args.true_z_dim else "random"

    rows = []
    for seed in range(args.seed_start, args.seed_start + args.n_sims):
        Z, W, X, G, beta = perturbvi.generate_sim(
            seed=seed, l_dim=args.true_l_dim, n_dim=args.n_dim, p_dim=args.p_dim,
            z_dim=args.true_z_dim, g_dim=args.g_dim, b_sparsity=args.b_sparsity,
        )
        G_sp = sparse.BCOO.fromdense(G)

        start = time.time()
        results = perturbvi.infer(
            X, G_sp, z_dim=args.z_dim, l_dim=args.l_dim,
            p_prior=args.p_prior, init=init, tau=args.tau,
            tol=args.tol, max_iter=args.max_iter, verbose=False,
        )
        running_time = time.time() - start

        beta_hat = results.params.mean_beta * results.params.p_hat.T
        lfsr = np.asarray(perturbvi.estimate_lfsr(results, draws=400, seed=0))
        overall = beta @ W

        beta_err = orthogonal(np.asarray(beta_hat), np.asarray(beta), scale=True, pad=True).error
        z_err = orthogonal(np.asarray(results.params.mean_z), np.asarray(Z), scale=True, pad=True).error
        w_err = orthogonal(np.asarray(results.W.T), np.asarray(W.T), scale=True, pad=True).error
        sen, spec = compute_sen_spec(overall, lfsr)

        rows.append({
            "beta_err": np.round(beta_err, 4),
            "z_err": np.round(z_err, 4),
            "w_err": np.round(w_err, 4),
            "sensitivity": np.round(sen, 4),
            "specificity": np.round(spec, 4),
            "running_time": running_time,
            "sim": seed,
        })
        print(f"seed {seed}: beta_err={beta_err:.4f} w_err={w_err:.4f} "
              f"sen={sen:.4f} ({running_time:.1f}s)")

    df = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / args.out
    df.to_csv(out_path, index=False)
    print(f"saved -> {out_path}")


if __name__ == "__main__":
    main()
