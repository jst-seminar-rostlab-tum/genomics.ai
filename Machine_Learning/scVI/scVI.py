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
    'vae_max_epochs': 20,
    'n_neighbors': 8,
    'max_epochs': 100,
    'unwanted_labels': ['leiden'],
    'model_name' : 'model.pt' # To see if the model already exists
}


def get_from_config(key):
    return config[key]



# python3.9 scVI.py --input data/pancreas_normalized.h5ad -t -q

def parse_args():
    parser = argparse.ArgumentParser(description='Run functions on scANVI using scArches')
    parser.add_argument('--input', help='.h5ad file containing the data to work on')
    parser.add_argument('--config', help='file containing all variable arguments needed to compute the result')
    parser.add_argument('-d', '--debug', help='print debug output', action='store_true')
    parser.add_argument('-t', '--train', help='train the scvi model', action='store_true')
    parser.add_argument('-q', '--query', help='execute query on reference dataset', action='store_true')
    #parser.add_argument('-p', '--predict', help='predict', action='store_true')
    return parser.parse_args()

def setup():        # Set up the warnings filter an the figure parameters
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)
    sc.settings.set_figure_params(dpi=200, frameon=False) # Set resolution/size, styling and format of figures.
    sc.set_figure_params(dpi=200, figsize=(4,4)) # https://scanpy.readthedocs.io/en/stable/generated/scanpy.set_figure_params.html
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7) # https://pytorch.org/docs/stable/generated/torch.set_printoptions.html

def pre_process_data(filename):         # Split the data in target and source, remove the sparsity
    condition_key = get_from_config('condition_key')
    cell_type_key = get_from_config('cell_type_key')
    target_conditions = get_from_config('target_conditions')

    # gdown.download('https://drive.google.com/uc?id=1ehxgfHTsMZXy6YzlFKGJOsBKQ5rrvMnd', 'pancreas.h5ad', quiet=False)

    adata_all = sc.read(filename)
    adata = adata_all.raw.to_adata()
    adata = remove_sparsity(adata)
    source_adata = adata[~adata.obs[condition_key].isin(target_conditions)].copy()
    target_adata = adata[adata.obs[condition_key].isin(target_conditions)].copy()
    return source_adata, target_adata

def create_scVI_model(source_adata, target_adata, logger, args): # if there is already a pretrained model, nothing happens otherwise a new one will be trained
    if path.exists((get_from_config('ref_path') + 'model.pt')):
        return sca.models.SCVI.load_query_data(
                    target_adata,
                    get_from_config('ref_path'),
                    freeze_dropout=True,
                )
    else:
        setup_anndata(source_adata)

        vae = get_model(source_adata)

        if args.debug:
            logger.debug(source_adata)
            logger.debug(target_adata)

        if args.train:
            vae.train()

        print(vae)

        compute_latent(vae, source_adata)

        vae.save(get_from_config('ref_path'), overwrite=True)

        return vae

def setup_anndata(anndata, filename):       # Just because it's prettier that way
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

def compute_query(anndata):     # trains the model on a query and saves the result
    model = sca.models.SCVI.load_query_data(
        anndata,
        get_from_config('ref_path'),
        freeze_dropout=True,
    )
    
    model.train(
        max_epochs=get_from_config('max_epochs'),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )
    query_latent = compute_latent(model, anndata)

    utils.save_umap_as_pdf(query_latent, 'figures/query.pdf', color=['batch', 'cell_type'])
    utils.write_latent_csv(query_latent, filename='/tmp/query.csv')
    return model, query_latent

def compute_full_latent(source_adata, target_adata, model):     # basically just takes to datasets, concatenates them and then computes the latent and saves the result
    full_latent = compute_latent(model, source_adata.concatenate(target_adata))
    utils.save_umap_as_pdf(full_latent, 'figures/both.pdf', color=['batch', 'cell_type'])
    utils.write_latent_csv(full_latent, filename='/tmp/both.csv')
    return full_latent

def main():
    args = parse_args()
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
    if args.config is None:
        logger.info('no config provided, using default config')
    filename = args.input
    show = True
    both = True

    setup()

    source_adata, target_adata = pre_process_data(filename)

    if args.debug:
        logger.debug(source_adata)
        logger.debug(target_adata)

    create_scVI_model(source_adata, target_adata, logger, args)

    model = None
    if args.query:
        model_query, query_latent = compute_query(target_adata)
        model = model_query
    
    full_latent = None

    if both:
        full_latent = compute_full_latent(source_adata, target_adata, model)



main()