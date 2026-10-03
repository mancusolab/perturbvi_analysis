from pathlib import Path

import anndata as ad
import jax.numpy as jnp
import pandas as pd
from jax import config
from jax.experimental import sparse

from perturbvi import PerturbData, estimate_lfsr, fit_screen, save_results

config.update("jax_enable_x64", True)
config.update("jax_default_matmul_precision", "highest")

Z_DIM, L_DIM, TAU, INIT = 20, 1000, 1, "pca"

DATA = Path("input")
MATRIX = DATA / "K562_essential_resid.h5ad"
GUIDE = DATA / "wide_df.csv"
BACKGROUND = DATA / "K562_essential_downstream_gene.tsv"
OUTPUT = Path("results")
DROP_COLS = ["non-targeting", "cell_barcode"]

adata = ad.read_h5ad(MATRIX)
guides = pd.read_csv(GUIDE, index_col=0)
guides = guides.drop(columns=DROP_COLS, errors="ignore")
background_genes = pd.read_csv(BACKGROUND, sep="\t")["gene_id"].tolist()

screen = PerturbData(
    X=jnp.asarray(adata.X, dtype=jnp.float64),
    G=sparse.bcoo_fromdense(jnp.asarray(guides.to_numpy(), dtype=jnp.float64)),
    gene_names=background_genes,
    perturbation_names=guides.columns.tolist(),
)

del adata, guides

fit = fit_screen(
    screen,
    z_dim=Z_DIM,
    l_dim=L_DIM,
    tau=TAU,
    init=INIT,
    p_prior=0.1,
    standardize=True,
    tol=1e-2,
    max_iter=1000,
)

save_results(fit, OUTPUT)
del fit

lfsr_bw = estimate_lfsr(OUTPUT)
lfsr_bw.to_csv(OUTPUT / "LFSR_BW.csv")
