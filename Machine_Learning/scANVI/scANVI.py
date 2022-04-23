import os
import warnings
from os import path
from scVI.scVI import create_scVI_model, set_config

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

config = None

def get_from_config(key):
    return config[key]

# url = 'https://drive.google.com/uc?id=1ehxgfHTsMZXy6YzlFKGJOsBKQ5rrvMnd'
# gdown.download(url, output, quiet=False)

def setup_modules():

    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)

    scanpy.settings.set_figure_params(dpi=200, frameon=False)
    scanpy.set_figure_params(dpi=200)
    scanpy.set_figure_params(figsize=(4, 4))

    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)


def pre_process_data(source_path, target_path):

    source_adata = scanpy.read(source_path)
    target_adata = scanpy.read(target_path)
    source_adata = remove_sparsity(source_adata)
    target_adata = remove_sparsity(target_adata)

    return source_adata, target_adata


def setup_anndata_for_scanvi(anndata):
    scarches.models.SCANVI.setup_anndata(anndata, unlabeled_category='Unknown',
                                         batch_key=get_from_config('condition_key'),
                                         labels_key=get_from_config('cell_type_key'))


def get_scanvi_from_scvi_model(scvi_model):
    return scarches.models.SCANVI.from_scvi_model(scvi_model, get_from_config('unlabeled_key'))


'''def get_scvi_model(anndata):
    return scarches.models.SCVI(
        anndata,
        n_layers=get_from_config('n_layers'),
        encode_covariates=get_from_config('encode_covariates'),
        deeply_inject_covariates=get_from_config('deeply_inject_covariates'),
        use_layer_norm=get_from_config('use_layer_norm'),
        use_batch_norm=get_from_config('use_batch_norm'),
    )'''


def get_latent(model, adata):

    reference_latent = scanpy.AnnData(model.get_latent_representation())
    reference_latent.obs["cell_type"] = adata.obs[get_from_config('cell_type_key')].tolist()
    reference_latent.obs["batch"] = adata.obs[get_from_config('condition_key')].tolist()

    scanpy.pp.neighbors(reference_latent, n_neighbors=get_from_config('n_neighbors'))
    scanpy.tl.leiden(reference_latent)
    scanpy.tl.umap(reference_latent)

    model.save(get_from_config('ref_path'), overwrite=True)

    return reference_latent


def predict(model, latent):
    latent.obs['predictions'] = model.predict()
    print("Acc: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
    return latent


def surgery(anndata, model_path, surgery_path, generate_output):

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
        max_epochs=100,
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )

    surgery_latent = get_latent(model, anndata)

    if generate_output:
        utils.save_umap_as_pdf(surgery_latent, 'figures/surgery.pdf', color=['batch', 'cell_type'])
        utils.write_latent_csv(surgery_latent, filename='/tmp/surgery.csv')

    model.save(surgery_path, overwrite=True)

    return model, surgery_latent


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

    figure = plt.figure(figsize=(8, 8), frameon=False)
    _ = plt.grid(False)
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


def compute_scANVI(configP, reference_dataset, query_dataset, model_path, surgery_path, generate_output):

    global config
    config = configP

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    both = False
    compare = False
    query = True
    predict = True

    setup_modules()

    source_adata, target_adata = pre_process_data(reference_dataset, query_dataset)

    set_config(configP) # sets the config in scVI

    vae = create_scVI_model(source_adata, target_adata, model_path) # kann man sehen ob es schon ein scVI oder ein scANVI model gibt?
    '''if args.debug:
        logger.debug(source_adata)
        logger.debug(target_adata)

    setup_anndata_for_scvi(source_adata)

    vae = get_scvi_model(source_adata)

    if args.debug:
        logger.debug(source_adata)
        logger.debug(target_adata)

    #if args.train:
    vae.train(max_epochs=get_from_config('scvi_max_epochs'))

    print(vae)'''

    # setup_anndata_for_scanvi(source_adata)

    scanvi = get_scanvi_from_scvi_model(vae)

    print("Labelled Indices: ", len(scanvi._labeled_indices))
    print("Unlabelled Indices: ", len(scanvi._unlabeled_indices))

    #print(scanvi)

    scanvi.train(max_epochs=get_from_config('scanvi_max_epochs'))
    reference_latent = get_latent(scanvi, source_adata)

    if generate_output:
        utils.save_umap_as_pdf(reference_latent, 'figures/reference.pdf', color=['batch', 'cell_type'])

    reference_latent = predict(scanvi, reference_latent)

    model = None
    if query:
        model_query, query_latent = query(target_adata)
        model = model_query

        if predict:
            predict_latent(predict(model_query, query_latent))

    #if args.surgery:
    model_surgery, surgery_latent = surgery(target_adata, model_path, surgery_path, generate_output)

    if predict:
        predict_latent(predict(model_surgery, surgery_latent))

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

