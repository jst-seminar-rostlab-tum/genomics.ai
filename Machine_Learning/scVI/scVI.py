import os
import os.path
from os import path

import warnings

import numpy as np
import gdown
import sys
import getopt
import torch
import tempfile
import utils
import logging, sys
import argparse
import scanpy as sc
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt

#config = None

def set_config(configP):
    global config
    config = configP

def get_from_config(key):
    
    return config[key]
    
# python3.9 scVI.py --input data/pancreas_normalized.h5ad -t -q


def setup():        # Set up the warnings filter an the figure parameters

    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)
    sc.settings.set_figure_params(dpi=200, frameon=False) # Set resolution/size, styling and format of figures.
    sc.set_figure_params(dpi=200, figsize=(4,4)) # https://scanpy.readthedocs.io/en/stable/generated/scanpy.set_figure_params.html
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7) # https://pytorch.org/docs/stable/generated/torch.set_printoptions.html


def pre_process_data(source_path, target_path):

    source_adata = sc.read(source_path)
    target_adata = sc.read(target_path)
    source_adata = remove_sparsity(source_adata)
    target_adata = remove_sparsity(target_adata)

    return source_adata, target_adata


def create_scVI_model(source_adata, target_adata, model_path): # if there is already a pretrained model, nothing happens otherwise a new one will be trained

    if path.exists(model_path):
        return sca.models.SCVI.load_query_data(
                    target_adata,
                    model_path,
                    freeze_dropout=True,
                )
    else:
        setup_anndata(source_adata)

        vae = get_model(source_adata)

        vae.train()

        compute_latent(vae, source_adata)

        vae.save(model_path, overwrite=True)

        return vae


def setup_anndata(anndata):       # Just because it's prettier that way
    sca.models.SCVI.setup_anndata(anndata, batch_key=get_from_config('condition_key'), labels_key=get_from_config('cell_type_key'))


def get_model(anndata):                     # Just because it's prettier that way
    return sca.models.SCVI(
        anndata,
        n_layers=get_from_config('n_layers'),
        encode_covariates=get_from_config('encode_covariates'),
        deeply_inject_covariates=get_from_config('deeply_inject_covariates'),
        use_layer_norm=get_from_config('use_layer_norm'),
        use_batch_norm=get_from_config('use_batch_norm'),
    )


def compute_latent(model, adata):       # computes the latent of a model with specific adata

    reference_latent = sc.AnnData(model.get_latent_representation(adata=adata))
    reference_latent.obs["cell_type"] = adata.obs[get_from_config('cell_type_key')].tolist()
    reference_latent.obs["batch"] = adata.obs[get_from_config('condition_key')].tolist()
    sc.pp.neighbors(reference_latent, n_neighbors=get_from_config('n_neighbors'))
    sc.tl.leiden(reference_latent)
    sc.tl.umap(reference_latent)
    sc.pl.umap(reference_latent,
           color=['batch', 'cell_type'],
           frameon=False,
           wspace=0.6,
           )
    return reference_latent


def compute_query(anndata, model_path, generate_output):     # trains the model on a query and saves the result

    model = sca.models.SCVI.load_query_data(
        anndata,
        model_path,
        freeze_dropout=True,
    )
    
    model.train(
        max_epochs=get_from_config('max_epochs'),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )
    query_latent = compute_latent(model, anndata)

    if generate_output:
        utils.save_umap_as_pdf(query_latent, 'data/figures/query.pdf', color=['batch', 'cell_type'])
        utils.write_latent_csv(query_latent, filename='/tmp/query.csv')

    return model


def compute_full_latent(source_adata, target_adata, model, generate_output):     # basically just takes to datasets, concatenates them and then computes the latent and saves the result

    full_latent = compute_latent(model, source_adata.concatenate(target_adata))

    if generate_output:
        utils.save_umap_as_pdf(full_latent, 'data/figures/both.pdf', color=['batch', 'cell_type'])
        utils.write_latent_csv(full_latent, filename='/tmp/both.csv')

    return full_latent

def compute_scVI(configP, reference_dataset, query_dataset, model_path, surgery_path, generate_output):

    global config
    config = configP

    setup()

    source_adata, target_adata = pre_process_data(reference_dataset, query_dataset)

    create_scVI_model(source_adata, target_adata, model_path)

    model = compute_query(target_adata, model_path, generate_output)

    compute_full_latent(source_adata, target_adata, model, generate_output)

    model.save(surgery_path, overwrite=True)


#main()