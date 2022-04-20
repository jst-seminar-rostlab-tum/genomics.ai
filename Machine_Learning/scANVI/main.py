import os
import warnings

import matplotlib
import scanpy
import scarches
from scarches.dataset.trvae.data_handling import remove_sparsity
from matplotlib import pyplot as plt
import numpy as np
import gdown
import sys
import getopt
import torch
import tempfile
import utils
import logging, sys
import argparse

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
    'unwanted_labels': ['leiden']
}


def get_from_config(key):
    return config[key]


# url = 'https://drive.google.com/uc?id=1ehxgfHTsMZXy6YzlFKGJOsBKQ5rrvMnd'
# gdown.download(url, output, quiet=False)


def parse_args():
    parser = argparse.ArgumentParser(description='Run functions on scANVI using scArches')
    parser.add_argument('--input', help='.h5ad file containing the data to work on')
    parser.add_argument('--config', help='file containing all variable arguments needed to compute the result')
    parser.add_argument('-d', '--debug', help='print debug output', action='store_true')
    parser.add_argument('-t', '--train', help='train the scvi model', action='store_true')
    parser.add_argument('-q', '--query', help='execute query on reference dataset', action='store_true')
    parser.add_argument('-p', '--predict', help='predict', action='store_true')
    return parser.parse_args()


def setup_modules():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)
    scanpy.settings.set_figure_params(dpi=200, frameon=False)
    scanpy.set_figure_params(dpi=200)
    scanpy.set_figure_params(figsize=(4, 4))
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)


def setup_anndata_for_scvi(anndata):
    scarches.models.SCVI.setup_anndata(anndata, batch_key=get_from_config('condition_key'), labels_key=get_from_config('cell_type_key'))


def setup_anndata_for_scanvi(anndata):
    scarches.models.SCANVI.setup_anndata(anndata, unlabeled_category='Unknown',
                                         batch_key=get_from_config('condition_key'), labels_key=get_from_config('cell_type_key'))


def get_scanvi_from_scvi_model(scvi_model):
    return scarches.models.SCANVI.from_scvi_model(scvi_model, get_from_config('unlabeled_key'), labels_key=get_from_config('cell_type_key'))


def get_scvi_model(anndata):
    return scarches.models.SCVI(
        anndata,
        n_layers=get_from_config('n_layers'),
        encode_covariates=get_from_config('encode_covariates'),
        deeply_inject_covariates=get_from_config('deeply_inject_covariates'),
        use_layer_norm=get_from_config('use_layer_norm'),
        use_batch_norm=get_from_config('use_batch_norm'),
    )


def get_latent(model, adata):
    reference_latent = scanpy.AnnData(model.get_latent_representation())
    reference_latent.obs["cell_type"] = adata.obs[get_from_config('cell_type_key')].tolist()
    reference_latent.obs["batch"] = adata.obs[get_from_config('condition_key')].tolist()
    scanpy.pp.neighbors(reference_latent, n_neighbors=get_from_config('n_neighbors'))
    scanpy.tl.leiden(reference_latent)
    scanpy.tl.umap(reference_latent)
    return reference_latent

def predict(model, latent):
    latent.obs['predictions'] = model.predict()
    print("Acc: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
    return latent

def query(anndata):
    model = scarches.models.SCANVI.load_query_data(
        anndata,
        get_from_config('ref_path'),
        freeze_dropout=True,
    )
    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []
    print("Labelled Indices: ", len(model._labeled_indices))
    print("Unlabelled Indices: ", len(model._unlabeled_indices))

    model.train(
        max_epochs=get_from_config('max_epochs'),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )
    query_latent = get_latent(model, anndata)

    utils.save_umap_as_pdf(query_latent, 'figures/query.pdf', color=['batch', 'cell_type'])
    utils.write_latent_csv(query_latent, filename='/tmp/query.csv')
    return model, query_latent


def predict_latent(latent):
    df = latent.obs.groupby(["cell_type", "predictions"]).size().unstack(fill_value=0)
    norm_df = df / df.sum(axis=0)

    figure = plt.figure(figsize=(8, 8), frameon=False, show=False)
    _ = plt.pcolor(norm_df)
    _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xlabel("Predicted")
    plt.ylabel("Observed")
    figure.savefig('figures/predict.png')


def both_adata(source_adata, target_adata):
    adata_full = source_adata.concatenate(target_adata)
    full_latent = scanpy.AnnData(scarches.models.SCANVI.get_latent_representation(adata=adata_full))
    full_latent.obs['cell_type'] = adata_full.obs[get_from_config('cell_type_key')].tolist()
    full_latent.obs['batch'] = adata_full.obs[get_from_config('condition_key')].tolist()
    scanpy.pp.neighbors(full_latent)
    scanpy.tl.leiden(full_latent)
    scanpy.tl.umap(full_latent)
    utils.save_umap_as_pdf(full_latent, 'figures/both.pdf', color=['batch', 'cell_type'])
    utils.write_latent_csv(full_latent, filename='/tmp/both.csv')
    return full_latent


def compare_adata(model, source_adata, target_adata, latent):
    adata_full = source_adata.concatenate(target_adata)
    latent.obs['predictions'] = model.predict(adata=adata_full)
    print("Acc: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
    scanpy.pp.neighbors(latent)
    scanpy.tl.leiden(latent)
    scanpy.tl.umap(latent)
    utils.save_umap_as_pdf(latent, 'figures/compare.pdf', color=["predictions", "cell_type"])
    utils.write_latent_csv(latent, filename='/tmp/compare.csv')


def main():
    args = parse_args()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
    if args.config is None:
        logger.info('no config provided, using default config')
    filename = args.input
    show = True
    both = False
    compare = False
    setup_modules()
    condition_key = get_from_config('condition_key')
    cell_type_key = get_from_config('cell_type_key')
    target_conditions = get_from_config('target_conditions')
    adata_all = scanpy.read(filename)
    adata = adata_all.raw.to_adata()
    adata = remove_sparsity(adata)
    source_adata = adata[~adata.obs[get_from_config('condition_key')].isin(target_conditions)].copy()
    target_adata = adata[adata.obs[get_from_config('condition_key')].isin(target_conditions)].copy()

    if args.debug:
        logger.debug(source_adata)
        logger.debug(target_adata)

    setup_anndata_for_scvi(source_adata)

    vae = get_scvi_model(source_adata)

    if args.debug:
        logger.debug(source_adata)
        logger.debug(target_adata)

    if args.train:
        vae.train()

    print(vae)

    #setup_anndata_for_scanvi(source_adata)

    scanvae = get_scanvi_from_scvi_model(vae)

    print("Labelled Indices: ", len(scanvae._labeled_indices))
    print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))

    print(scanvae)

    scanvae.train(max_epochs=get_from_config('scanvae_max_epochs'))
    reference_latent = get_latent(scanvae, source_adata)

    if show:
        utils.save_umap_as_pdf(reference_latent, 'figures/reference.pdf', color=['batch', 'cell_type'])

    reference_latent = predict(scanvae, reference_latent)

    model = None
    if args.query:
        model_query, query_latent = query(target_adata)
        model = model_query

        if args.predict:
            predict_latent(query_latent)

    full_latent = None

    if both:
        full_latent = both_adata(source_adata, target_adata)

    if compare:
        if full_latent is None:
            full_latent = both_adata(source_adata, target_adata)
        if model is None:
            model_query, query_latent = query(target_adata)
            model = model_query
        compare_adata(model, source_adata, target_adata, full_latent)


if __name__ == "__main__":
    main()
