#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from time import time
import argparse as ap
import sys

import perturbvi
import numpy as np
import pandas as pd
import pickle
from perturbvi import utils

import jax.numpy as jnp
from jax import config
from jax.experimental import sparse
import jax.random as rdm

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from adjustText import adjust_text
import tools
import csv 
import psutil
import os

config.update("jax_enable_x64", True)
config.update("jax_default_matmul_precision", "highest")


path = "/Users/camellia/Project/perturbvi/data"
data = pd.read_csv(f"{path}/luhmes_exp.csv",index_col=0)
G = pd.read_csv(f"{path}/luhmes_G.csv",index_col=0)
G_reduce = G.drop(columns=["Nontargeting"])
g_sp = sparse.bcoo_fromdense(G)
g_reduce_sp = sparse.bcoo_fromdense(G_reduce)

X = jnp.asarray(data)
print(f"Begin inference on local cpu using raw data")
start = time()
results = perturbvi.infer(X,
                z_dim=12,
                l_dim=400,
                G = g_reduce_sp,
                A = None,
                p_prior = 0.5,
                standardize = False,
                init="random",
                tau = 10,
                max_iter = 500,
                tol=1e-2)
end = time()
print(f"inference Finished in {end-start}")
perturbvi.utils.pip_analysis(results.pip,rho=0.9,rho_prime=0.10)

# test memory use
def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem_before = process.memory_info().rss / 1024 / 1024

    # Function call
    results = perturbvi.infer(X,
                       z_dim=15,
                       l_dim=450,
                       G=g_reduce_sp,
                       A=None,
                       p_prior=0.5,
                       standardize=True,
                       tau=800,
                       max_iter=400,
                       tol=1e-2)

    mem_after = process.memory_info().rss / 1024 / 1024
    print(f"Memory usage: {mem_after - mem_before:.2f} MB")

memory_usage_psutil()




####################################

params = results.params
pip = results.pip
np.sum(results.pve)
#SuSiE PCA also show 4 factor with 0 effects
# results_sp = sp.infer.susie_pca(X,
#                 z_dim=10,
#                 l_dim=300,
#                 A = None,
#                 standardize = True,
#                 tau = 10,
#                 max_iter = 400,
#                 tol=1e-2)

# sp.utils.pip_analysis(results_sp.pip,rho=0.9,rho_prime=0.05)

#save results
#sp.io.save_results(results,path = "/Users/dongyuan/Documents/Project/gSuSiEPCA/luhmes/luhmes_results/Perturb_vi")

# first set K=20 and L=300, last 5 factors does not contribute; pve = 1.99%
# then set K = 15 and L=400, Last 2 factors does not contribute; pve = 2.11%
# set K = 12 and L = 400; current version; 166 seconds

#load gene symbol
with open(f"{path}/luhmes_gene_symbol.csv", mode='r', encoding='utf-8') as file:
    reader = csv.reader(file)
    gene_symbol = [row[0] for row in reader]
#column name
z_dim, p_dim = params.W.shape
g_dim, z_dim = params.mean_beta.shape
n_dim, z_dim = params.mean_z.shape
column_names_w = ['w' + str(i) for i in range(z_dim)]
column_names_pip = ['pip' + str(i) for i in range(z_dim)]
column_names_z = ['z' + str(i) for i in range(z_dim)]
column_names_b = ['b' + str(i) for i in range(z_dim)]
pip_df = pd.DataFrame(pip.T, columns=column_names_w, index=gene_symbol)
perturb_gene_list = G_reduce.columns.tolist()

path_results = "/Users/camellia/Project/perturbvi_analysis/results/luhmes"
pip_df.to_csv(f"{path_results}/pip_df.csv")

perturb_degs = tools.find_top_genes(pip_df,0.95)

np.asarray(perturb_degs["w0"])
np.asarray(perturb_degs["w1"])

np.asarray(perturb_degs["w3"])
np.asarray(perturb_degs["w4"])
np.asarray(perturb_degs["w5"])

np.asarray(perturb_degs["w7"])


beta_sparse = params.mean_beta * params.p_hat.T
beta_sparse_df = pd.DataFrame(beta_sparse, columns=column_names_b, index=G_reduce.columns.tolist())
beta_sparse_df.to_csv(f"{path_results}/beta_target.csv")
div = sns.diverging_palette(250, 10, as_cmap=True, center = "light")
div_prob = sns.light_palette("navy", n_colors=256, reverse=False, as_cmap=True)
sns.heatmap(beta_sparse_df ,cmap=div,center=0,vmin=-np.max(np.abs(beta_sparse_df)),vmax=np.max(np.abs(beta_sparse_df)))
plt.xlabel("Factors")
plt.ylabel("Perturbations")

np.asarray(perturb_degs["w2"])
np.asarray(perturb_degs["w6"])
np.asarray(perturb_degs["w10"])

# compute lfsr
lfsr = utils.compute_lfsr(params)
# subset if rows contain a value that is < 0.05
lfsr_df = pd.DataFrame(lfsr, index=G_reduce.columns.tolist())
lfsr_df.to_csv(f"{path_results}/lfsr_df.csv")
sig_df = lfsr_df[(lfsr_df < 0.05).sum(axis=1)>1]
# test the new lfsr function
lfsr_key = rdm.PRNGKey(0)
lfsr_new = utils.compute_lfsr(lfsr_key, params)
lfsr_new.block_until_ready()
lfsr_new_np = np.array(lfsr_new)  # Convert JAX array to NumPy
lfsr_new_df = pd.DataFrame(lfsr_new_np, index=G_reduce.columns.tolist())
#lfsr_new_df.to_csv(f"{path_results}/lfsr_new_df.csv")

