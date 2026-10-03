"""Simulate missing targets: drop a fraction of guide columns and fit PerturbVI
with both the reduced and the correct design.

Scenario (80 seeds, seeds 0..79; n=3000, p=4000, K=4, L=150, g=100, b=0.20):
  remove_prop: 0.1, 0.2, 0.3, 0.4, 0.5

Usage:
  python simulation_mis_g.py --remove_prop 0.1 --out perturbvi_misg_10.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from jax import config, random
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
        description="Simulate missing targets (dropped guides)."
    )
    p.add_argument("--remove_prop", type=float, default=0.4,
                   help="fraction of perturbation columns to drop")
    p.add_argument("--n_sims", type=int, default=80)
    p.add_argument("--seed_start", type=int, default=0)
    p.add_argument("--out", type=str, default="simulation_mis_g.csv")
    return p.parse_args()


def main():
    args = parse_args()

    # fixed scenario constants matching the archived mis_g runs
    n_dim = 3000
    p_dim = 4000
    real_l_dim = 150
    l_dim = 150
    z_dim = 4
    g_dim = 100
    b_sparsity = 0.20

    rows = []
    for seed in range(args.seed_start, args.seed_start + args.n_sims):
        Z, W, X, G, beta = perturbvi.generate_sim(
            seed=seed, l_dim=real_l_dim, n_dim=n_dim, p_dim=p_dim,
            z_dim=z_dim, g_dim=g_dim, b_sparsity=b_sparsity,
        )

        # Drop a fixed fraction of perturbation columns, deterministic per seed.
        remove_idx = np.asarray(
            random.choice(random.PRNGKey(seed), g_dim,
                          shape=(int(g_dim * args.remove_prop),), replace=False)
        )
        G_sub = np.delete(G, remove_idx, axis=1)
        beta_sub = np.delete(beta, remove_idx, axis=0)

        G_sp = sparse.BCOO.fromdense(G)
        G_sub_sp = sparse.BCOO.fromdense(G_sub)

        results = perturbvi.infer(
            X, G_sub_sp, z_dim=z_dim, l_dim=l_dim,
            p_prior=0.5, init="random", tau=1, tol=1e-3, max_iter=500,
            verbose=False,
        )
        results_correct = perturbvi.infer(
            X, G_sp, z_dim=z_dim, l_dim=l_dim,
            p_prior=0.5, init="random", tau=1, tol=1e-3, max_iter=500,
            verbose=False,
        )

        beta_hat = results.params.mean_beta * results.params.p_hat.T
        lfsr = np.asarray(perturbvi.estimate_lfsr(results, draws=400, seed=0))
        lfsr_correct = np.asarray(
            perturbvi.estimate_lfsr(results_correct, draws=400, seed=0)
        )
        lfsr_correct_sub = np.delete(lfsr_correct, remove_idx, axis=0)

        overall_sub = beta_sub @ W

        beta_err = orthogonal(np.asarray(beta_hat), np.asarray(beta_sub), scale=True, pad=True).error
        z_err = orthogonal(np.asarray(results.params.mean_z), np.asarray(Z), scale=True, pad=True).error
        w_err = orthogonal(np.asarray(results.W.T), np.asarray(W.T), scale=True, pad=True).error
        sen, spec = compute_sen_spec(overall_sub, lfsr)
        sen_correct, spec_correct = compute_sen_spec(overall_sub, lfsr_correct_sub)

        rows.append({
            "beta_err": np.round(beta_err, 4),
            "z_err": np.round(z_err, 4),
            "w_err": np.round(w_err, 4),
            "sensitivity": np.round(sen, 4),
            "specificity": np.round(spec, 4),
            "sensitivity_correct": np.round(sen_correct, 4),
            "specificity_correct": np.round(spec_correct, 4),
            "sim": seed,
        })
        print(f"seed {seed}: sen={sen:.4f} sen_correct={sen_correct:.4f}")

    df = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / args.out
    df.to_csv(out_path, index=False)
    print(f"saved -> {out_path}")


if __name__ == "__main__":
    main()
