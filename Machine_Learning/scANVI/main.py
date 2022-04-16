import os
import warnings
import scanpy as sc
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import sys
import getopt
import torch

config = {
    'condition_key': 'study',
    'cell_type_key': 'cell_type',
    'target_conditions': ['Pancreas CelSeq2', 'Pancreas SS2'],
    'ref_path': './ref_model/',
    'n_layers': 2,
    'encode_covariates': True,
    'deeply_inject_covariates': False,
    'use_layer_norm': 'both',
    'use_batch_norm': 'none',
    'unlabeled_key': 'Unknown',
    'scanvae_max_epochs': 20,
    'n_neighbors': 8,
    'max_epochs': 100,
}

def get_from_config(key):
    return config[key]

# url = 'https://drive.google.com/uc?id=1ehxgfHTsMZXy6YzlFKGJOsBKQ5rrvMnd'
# gdown.download(url, output, quiet=False)

def parse_arguments():
    input = ''
    config = ''
    try:
        opts, args = getopt.getopt(sys.argv, "hi:c:", ["input=", "config="])
    except getopt.GetoptError:
        print(__file__ + ' -i <inputfile> -c <configfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(__file__ + ' -i <inputfile> -c <configfile>')
            sys.exit()
        elif opt in ("-i", "--input"):
            input = arg
        elif opt in ("-o", "--config"):
            config = arg
    return input, config


def main():
    input, config = parse_arguments()
    os.chdir('../')
    if config is None:
        print('no config provided, using default config')
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)
    filename = 'pancreas.h5ad'
    show = True
    train = True
    query = True
    file_suffix = 6
    sc.settings.set_figure_params(dpi=200, frameon=False)
    sc.set_figure_params(dpi=200)
    sc.set_figure_params(figsize=(4, 4))
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
    condition_key = get_from_config('condition_key')
    cell_type_key = get_from_config('cell_type_key')
    target_conditions = get_from_config('target_conditions')
    ref_path = get_from_config('ref_path')
    adata_all = sc.read(filename)
    adata = adata_all.raw.to_adata()
    adata = remove_sparsity(adata)
    source_adata = adata[~adata.obs[condition_key].isin(target_conditions)].copy()
    target_adata = adata[adata.obs[condition_key].isin(target_conditions)].copy()
    sca.models.SCVI.setup_anndata(source_adata, batch_key=condition_key, labels_key=cell_type_key)
    vae = sca.models.SCVI(
        source_adata,
        n_layers=get_from_config('n_layers'),
        encode_covariates=get_from_config('encode_covariates'),
        deeply_inject_covariates=get_from_config('deeply_inject_covariates'),
        use_layer_norm=get_from_config('use_layer_norm'),
        use_batch_norm=get_from_config('use_batch_norm'),
    )
    if train:
        vae.train()
    scanvae = sca.models.SCANVI.from_scvi_model(vae, get_from_config('unlabeled_key'))
    #sca.models.SCANVI.setup_anndata(source_adata, batch_key=condition_key, labels_key=cell_type_key, unlabeled_category=get_from_config('unlabeled_key'))
    #scanvae = sca.models.SCANVI(source_adata)
    #scanvae = sca.models.SCANVI(source_adata)
    print("Labelled Indices: ", len(scanvae._labeled_indices))
    print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))
    scanvae.train(max_epochs=get_from_config('scanvae_max_epochs'))
    reference_latent = sc.AnnData(scanvae.get_latent_representation())
    reference_latent.obs["cell_type"] = source_adata.obs[cell_type_key].tolist()
    reference_latent.obs["batch"] = source_adata.obs[condition_key].tolist()
    sc.pp.neighbors(reference_latent, n_neighbors=get_from_config('n_neighbors'))
    sc.tl.leiden(reference_latent)
    sc.tl.umap(reference_latent)

    if show:
        sc.pl.umap(reference_latent,
                   color=['batch', 'cell_type'],
                   frameon=False,
                   wspace=0.6,
                   )

    reference_latent.obs['predictions'] = scanvae.predict()
    print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))
    if query:
        model = sca.models.SCANVI.load_query_data(
            target_adata,
            ref_path,
            freeze_dropout=True,
        )
        model._unlabeled_indices = np.arange(target_adata.n_obs)
        model._labeled_indices = []
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", len(model._unlabeled_indices))

        model.train(
            max_epochs=get_from_config('max_epochs'),
            plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10,
        )
        query_latent = sc.AnnData(model.get_latent_representation())
        query_latent.obs['cell_type'] = target_adata.obs[cell_type_key].tolist()
        query_latent.obs['batch'] = target_adata.obs[condition_key].tolist()
        sc.pp.neighbors(query_latent)
        sc.tl.leiden(query_latent)
        sc.tl.umap(query_latent)
        plt.figure()
        if show:
            sc.pl.umap(
                query_latent,
                color=["batch", "cell_type"],
                frameon=False,
                wspace=0.6,
            )
        final = query_latent.obs.drop(columns=['leiden'])
        final["x"] = list(map(lambda p: p[0], query_latent.obsm["X_umap"]))
        final["y"] = list(map(lambda p: p[1], query_latent.obsm["X_umap"]))
        final.to_csv('result_query' + str(file_suffix) + '.csv')
    final = reference_latent.obs.drop(columns=['leiden'])
    final["x"] = list(map(lambda p: p[0], reference_latent.obsm["X_umap"]))
    final["y"] = list(map(lambda p: p[1], reference_latent.obsm["X_umap"]))
    final.to_csv('result_reference' + str(file_suffix) + '.csv')


if __name__ == "__main__":
    main()
