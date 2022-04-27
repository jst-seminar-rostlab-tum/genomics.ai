import os
import warnings
import scanpy as sc
import anndata
import torch
import scarches as sca
import matplotlib.pyplot as plt
import numpy as np
import scvi as scv
import pandas as pd
from os.path import exists

def main():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)

    sc.settings.set_figure_params(dpi=200, frameon=False)
    sc.set_figure_params(dpi=200)
    sc.set_figure_params(figsize=(4, 4))
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

    adata_ref = scv.data.pbmcs_10x_cite_seq()

    adata_query = scv.data.dataset_10x("pbmc_10k_v3")
    adata_query.obs["batch"] = "PBMC 10k (RNA only)"
    # put matrix of zeros for protein expression (considered missing)
    pro_exp = adata_ref.obsm["protein_expression"]
    data = np.zeros((adata_query.n_obs, pro_exp.shape[1]))
    adata_query.obsm["protein_expression"] = pd.DataFrame(columns=pro_exp.columns, index=adata_query.obs_names, data = data)

    adata_full = anndata.concat([adata_ref, adata_query])

    adata_ref = adata_full[np.logical_or(adata_full.obs.batch == "PBMC5k", adata_full.obs.batch == "PBMC10k")].copy()
    adata_query = adata_full[adata_full.obs.batch == "PBMC 10k (RNA only)"].copy()

    sc.pp.highly_variable_genes(
        adata_ref,
        n_top_genes=4000,
        flavor="seurat_v3",
        batch_key="batch",
        subset=True,
    )

    adata_query = adata_query[:, adata_ref.var_names].copy()

    sca.models.TOTALVI.setup_anndata(
        adata_ref,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression"
    )

    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
    )

    dir_path = 'saved_model/'
    file_exists = exists(dir_path + "model.pt")

    vae_ref = None
    if file_exists:
        vae_ref = sca.models.TOTALVI.load(adata=adata_ref, dir_path=dir_path)
    else:
        vae_ref = sca.models.TOTALVI(
            adata_ref,
            **arches_params
        )
        vae_ref.train()

    adata_ref.obsm["X_totalVI"] = vae_ref.get_latent_representation()
    sc.pp.neighbors(adata_ref, use_rep="X_totalVI")
    sc.tl.umap(adata_ref, min_dist=0.4)

    sc.pl.umap(
        adata_ref,
        color=["batch"],
        frameon=False,
        ncols=1,
        title="Reference",
        save=True
    )

    os.rename('figures/umap.pdf', "figures/firstumap.pdf")

    dir_path = "saved_model/"
    vae_ref.save(dir_path, overwrite=True)

    # ...


    dir_path_2 = "saved_model/" # change to saved_model_2 if different models are wanted
    vae_q_file = exists(dir_path_2 + "model.pt")
    vae_q = None
    #if vae_q_file:
    #    vae_q = sca.models.TOTALVI.load(adata=adata_query, dir_path=dir_path) #, extend_categories=True
    #else:
    vae_q = sca.models.TOTALVI.load_query_data(
        adata_query,
        dir_path,
        freeze_expression=True
    )
    vae_q.train(200, plan_kwargs=dict(weight_decay=0.0), use_gpu=True)

    adata_query.obsm["X_totalVI"] = vae_q.get_latent_representation()
    sc.pp.neighbors(adata_query, use_rep="X_totalVI")
    sc.tl.umap(adata_query, min_dist=0.4)

    vae_q.save(dir_path_2, overwrite=True)


    _, imputed_proteins = vae_q.get_normalized_expression(
        adata_query,
        n_samples=25,
        return_mean=True,
        transform_batch=["PBMC10k", "PBMC5k"],
    )

    adata_query.obs = pd.concat([adata_query.obs, imputed_proteins], axis=1)

    sc.pl.umap(
        adata_query,
        color=imputed_proteins.columns,
        frameon=False,
        ncols=3,
        save=True
    )

    os.rename('figures/umap.pdf', "figures/secondumap.pdf")


    adata_full_new = adata_query.concatenate(adata_ref, batch_key="none")
    adata_full_new.obsm["X_totalVI"] = vae_q.get_latent_representation(adata_full_new)
    sc.pp.neighbors(adata_full_new, use_rep="X_totalVI")
    sc.tl.umap(adata_full_new, min_dist=0.3)

    _, imputed_proteins_all = vae_q.get_normalized_expression(
        adata_full_new,
        n_samples=25,
        return_mean=True,
        transform_batch=["PBMC10k", "PBMC5k"],
    )

    for i, p in enumerate(imputed_proteins_all.columns):
        adata_full_new.obs[p] = imputed_proteins_all[p].to_numpy().copy()


    perm_inds = np.random.permutation(np.arange(adata_full_new.n_obs))
    sc.pl.umap(
        adata_full_new[perm_inds],
        color=["batch"],
        frameon=False,
        ncols=1,
        title="Reference and query",
        save=True
    )

    os.rename('figures/umap.pdf', "figures/thirdumap.pdf")

    ax = sc.pl.umap(
        adata_full_new,
        color="batch",
        groups=["PBMC 10k (RNA only)"],
        frameon=False,
        ncols=1,
        title="Reference and query",
        alpha=0.4,
        save=True
    )
    os.rename('figures/umap.pdf', "figures/forthumap.pdf")

    sc.pl.umap(
        adata_full_new,
        color=imputed_proteins_all.columns,
        frameon=False,
        ncols=3,
        vmax="p99",
        save=True
    )
    os.rename('figures/umap.pdf', "figures/fifthumap.pdf")

if __name__ == "__main__":
    main()
