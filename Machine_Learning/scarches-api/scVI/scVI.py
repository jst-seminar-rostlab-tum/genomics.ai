import os
import warnings

import scanpy as sc
import scarches as sca
import scvi
import torch
from scarches.dataset.trvae.data_handling import remove_sparsity

from utils import utils, parameters
import sys
import tempfile
import scvi


# def utils.get_from_config(configuration, key):
#     if key in configuration:
#         return configuration[key]
#     return None


# python3.9 scVI.py --input data/pancreas_normalized.h5ad -t -q


def setup():
    """
    Set up the warnings filter and the figure parameters
    :return:
    """
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)
    #  Set resolution/size, styling and format of figures.
    sc.settings.set_figure_params(dpi=200, frameon=False)
    # https://scanpy.readthedocs.io/en/stable/generated/scanpy.set_figure_params.html
    sc.set_figure_params(dpi=200, figsize=(4, 4))
    # https://pytorch.org/docs/stable/generated/torch.set_printoptions.html
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)


def get_pretrained_scVI_model(anndata, configuration):
    """
    returns pretrained and saved scvi model
    :param anndata: query data to be used on the model
    :param configuration: configuration containing the name of the atlas
    :return: scarches SCVI model
    """
    return sca.models.SCVI.load_query_data(
        anndata,
        'assets/scVI/' + str(utils.get_from_config(configuration, parameters.ATLAS)) + '/',
        freeze_dropout=True,
    )


def create_scVI_model(source_adata, target_adata, configuration):
    """
    if there is already a pretrained model, nothing happens otherwise a new one will be trained
    :param source_adata: reference data
    :param target_adata: query data
    :return:
    """
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        print('use_pretrained is ' + str(utils.get_from_config(configuration, parameters.USE_PRETRAINED_SCVI_MODEL)),
              file=sys.stderr)
    if utils.get_from_config(configuration, parameters.USE_PRETRAINED_SCVI_MODEL):
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            print('use pretrained scvi model', file=sys.stderr)
        # os.mkdir('scvi_model')
        # utils.fetch_file_from_s3(utils.get_from_config(configuration, parameters.PRETRAINED_MODEL_PATH), 'assets/scVI/model.pt')

        
        return get_pretrained_scVI_model(target_adata, configuration), None
    else:
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            print('do not use pretrained scvi model', file=sys.stderr)
        setup_anndata(source_adata, configuration)
        vae = get_model(source_adata, configuration)
        vae.train(max_epochs=utils.get_from_config(configuration, parameters.SCVI_MAX_EPOCHS),
                  use_gpu=utils.get_from_config(configuration, parameters.USE_GPU))
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.write_adata_to_csv(vae, source_adata, key='scvi-source-adata-post-first-training.csv')
            except Exception as e:
                print(e, file=sys.stderr)
        reference_latent = compute_latent(vae, source_adata, configuration)
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.write_adata_to_csv(vae, source_adata, key='scvi-reference-latent-post-first-training.csv')
            except Exception as e:
                print(e, file=sys.stderr)
        tempdir = tempfile.mkdtemp()
        vae.save(tempdir, overwrite=True)
        
        print(os.listdir(tempdir), file=sys.stderr)
        # utils.store_file_in_s3(tempdir + '/model.pt', utils.get_from_config(configuration, parameters.RESULTING_MODEL_PATH))
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.store_file_in_s3(tempdir + '/model.pt', 'scvi-model-after-first-training.pt')
            except Exception as e:
                print(e, file=sys.stderr)
        utils.delete_file(tempdir + '/model.pt')
        os.removedirs(tempdir)
        return vae, reference_latent


def setup_anndata(anndata, configuration):
    """
    wrapper around setup_anndata
    :param anndata:
    """
    sca.models.SCVI.setup_anndata(anndata, batch_key=utils.get_from_config(configuration, parameters.CONDITION_KEY),
                                  labels_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY))


def get_model(anndata, configuration):
    """
    wrapper around creating a SCVI model using the given configuration
    :param anndata:
    :return:
    """
    return sca.models.SCVI(
        anndata,
        n_layers=utils.get_from_config(configuration, parameters.NUMBER_OF_LAYERS),
        encode_covariates=utils.get_from_config(configuration, parameters.ENCODE_COVARIATES),
        deeply_inject_covariates=utils.get_from_config(configuration, parameters.DEEPLY_INJECT_COVARIATES),
        use_layer_norm=utils.get_from_config(configuration, parameters.USE_LAYER_NORM),
        use_batch_norm=utils.get_from_config(configuration, parameters.USE_BATCH_NORM),
    )


def compute_latent(model, adata, configuration):
    """
    computes the latent of a model with specific adata
    :param model:
    :param adata:
    :return:
    """
    reference_latent = sc.AnnData(model.get_latent_representation(adata=adata))
    reference_latent.obs[utils.get_from_config(configuration, parameters.CELL_TYPE_KEY)] = adata.obs[
        utils.get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
    reference_latent.obs[utils.get_from_config(configuration, parameters.CONDITION_KEY)] = adata.obs[
        utils.get_from_config(configuration, parameters.CONDITION_KEY)].tolist()
    sc.pp.neighbors(reference_latent, n_neighbors=utils.get_from_config(configuration, parameters.NUMBER_OF_NEIGHBORS))
    sc.tl.leiden(reference_latent)
    sc.tl.umap(reference_latent)

    return reference_latent


def compute_query(pretrained_model, anndata, reference_latent, source_adata, configuration):
    """
    trains the model on a query and saves the result
    :param anndata:
    :return:
    """
    model = sca.models.SCVI.load_query_data(
        anndata,
        'assets/scVI/' + str(utils.get_from_config(configuration, parameters.ATLAS)) + '/',
        freeze_dropout=True,
    )
    if utils.get_from_config(configuration, parameters.ATLAS) == 'human_lung':
        surgery_epochs = 500
        train_kwargs_surgery = {
            "early_stopping": True,
            "early_stopping_monitor": "elbo_train",
            "early_stopping_patience": 10,
            "early_stopping_min_delta": 0.001,
            "plan_kwargs": {"weight_decay": 0.0},
        }
        model.train(
            max_epochs=surgery_epochs,
            **train_kwargs_surgery,
            use_gpu=utils.get_from_config(configuration, parameters.USE_GPU)
        )
    else:
        model.train(
            max_epochs=utils.get_from_config(configuration, parameters.SCVI_QUERY_MAX_EPOCHS),
            plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10,
            use_gpu=utils.get_from_config(configuration, parameters.USE_GPU)
        )
    tempdir = tempfile.mkdtemp()
    model.save(tempdir, overwrite=True)
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'scvi-model-after-query-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
    utils.delete_file(tempdir + '/model.pt')
    os.removedirs(tempdir)
    
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            if reference_latent is not None:
                utils.write_latent_csv(reference_latent, key='reference-latent-post-query-training.csv')
            utils.write_adata_to_csv(model, anndata, key='query-adata-post-query-training.csv',
                                     cell_type_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                     condition_key=utils.get_from_config(configuration, parameters.CONDITION_KEY))
            utils.write_adata_to_csv(model, source_adata, key='source-adata-post-query-training.csv',
                                     cell_type_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                     condition_key=utils.get_from_config(configuration, parameters.CONDITION_KEY))
        except Exception as e:
            print(e, file=sys.stderr)

    query_latent = compute_latent(model, anndata, configuration)
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_latent_csv(query_latent, key='query-latent-post-query-training.csv')
            if reference_latent is not None:
                utils.write_combined_csv(reference_latent, query_latent, key='combined-latents-after-query.csv')
        except Exception as e:
            print(e, file=sys.stderr)
    if utils.get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(query_latent, 'data/figures/query.pdf', color=['batch', 'cell_type'])

    utils.write_full_adata_to_csv(model, source_adata, anndata,
                                  key=utils.get_from_config(configuration, parameters.OUTPUT_PATH),
                                  cell_type_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                  condition_key=utils.get_from_config(configuration, parameters.CONDITION_KEY), configuration=configuration)

    return model


def compute_scVI(configuration):
    setup()
    source_adata, target_adata = utils.pre_process_data(configuration)
    print(source_adata)
    print(target_adata)
    model, reference_latent = create_scVI_model(source_adata, target_adata, configuration)
    model = compute_query(model, target_adata, reference_latent, source_adata, configuration)
    # Saving of the pre-trained models on an organization level follows below

    # compute_full_latent(source_adata, target_adata, model)
    # model.save('resulting_model', overwrite=True)
    # utils.store_file_in_s3('resulting_model/model.pt', utils.get_from_config(configuration, parameters.RESULTING_MODEL_PATH) + '_new')
