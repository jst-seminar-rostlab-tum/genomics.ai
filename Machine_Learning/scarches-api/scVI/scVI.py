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


def get_from_config(configuration, key):
    if key in configuration:
        return configuration[key]
    return None


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


def pre_process_data(configuration):
    source_adata = utils.read_h5ad_file_from_s3(get_from_config(configuration, parameters.REFERENCE_DATA_PATH))
    target_adata = utils.read_h5ad_file_from_s3(get_from_config(configuration, parameters.QUERY_DATA_PATH))
    source_adata = remove_sparsity(source_adata)
    target_adata = remove_sparsity(target_adata)

    return source_adata, target_adata


def get_pretrained_scVI_model(anndata):
    return sca.models.SCVI.load_query_data(
        anndata,
        'assets/scVI/',
        freeze_dropout=True,
    )


def create_scVI_model(source_adata, target_adata, configuration):
    """
    if there is already a pretrained model, nothing happens otherwise a new one will be trained
    :param source_adata:
    :param target_adata:
    :return:
    """
    if get_from_config(configuration, parameters.DEV_DEBUG):
        print('use_pretrained is ' + str(get_from_config(configuration, parameters.USE_PRETRAINED_SCVI_MODEL)),
              file=sys.stderr)
    if get_from_config(configuration, parameters.USE_PRETRAINED_SCVI_MODEL):
        if get_from_config(configuration, parameters.DEV_DEBUG):
            print('use pretrained scvi model', file=sys.stderr)
        # os.mkdir('scvi_model')
        # utils.fetch_file_from_s3(get_from_config(configuration, parameters.PRETRAINED_MODEL_PATH), 'assets/scVI/model.pt')
        return get_pretrained_scVI_model(target_adata), None
    else:
        if get_from_config(configuration, parameters.DEV_DEBUG):
            print('do not use pretrained scvi model', file=sys.stderr)
        setup_anndata(source_adata, configuration)
        vae = get_model(source_adata, configuration)
        vae.train(max_epochs=get_from_config(configuration, parameters.SCVI_MAX_EPOCHS),
                  use_gpu=get_from_config(configuration, parameters.USE_GPU))

        #TODO check predict cell


        if get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.write_adata_to_csv(vae, source_adata, key='scvi-source-adata-post-first-training.csv')
            except Exception as e:
                print(e, file=sys.stderr)
        reference_latent = compute_latent(vae, source_adata, configuration)
        if get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.write_adata_to_csv(vae, source_adata, key='scvi-reference-latent-post-first-training.csv')
            except Exception as e:
                print(e, file=sys.stderr)
        tempdir = tempfile.mkdtemp()
        vae.save(tempdir, overwrite=True)
        # TODO check if we need to to this
        print(os.listdir(tempdir), file=sys.stderr)
        # utils.store_file_in_s3(tempdir + '/model.pt', get_from_config(configuration, parameters.RESULTING_MODEL_PATH))
        if get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.store_file_in_s3(tempdir + '/model.pt', 'scvi-model-after-first-training.pt')
            except Exception as e:
                print(e, file=sys.stderr)
        utils.delete_file(tempdir + '/model.pt')
        os.removedirs(tempdir)
        return vae, reference_latent


def setup_anndata(anndata, configuration):
    """
    Just because it's prettier that way
    :param anndata:
    :return:
    """
    sca.models.SCVI.setup_anndata(anndata, batch_key=get_from_config(configuration, parameters.CONDITION_KEY),
                                  labels_key=get_from_config(configuration, parameters.CELL_TYPE_KEY))


def get_model(anndata, configuration):
    """
    Just because it's prettier that way
    :param anndata:
    :return:
    """
    return sca.models.SCVI(
        anndata,
        n_layers=get_from_config(configuration, parameters.NUMBER_OF_LAYERS),
        encode_covariates=get_from_config(configuration, parameters.ENCODE_COVARIATES),
        deeply_inject_covariates=get_from_config(configuration, parameters.DEEPLY_INJECT_COVARIATES),
        use_layer_norm=get_from_config(configuration, parameters.USE_LAYER_NORM),
        use_batch_norm=get_from_config(configuration, parameters.USE_BATCH_NORM),
    )


def compute_latent(model, adata, configuration):
    """
    computes the latent of a model with specific adata
    :param model:
    :param adata:
    :return:
    """
    reference_latent = sc.AnnData(model.get_latent_representation(adata=adata))
    reference_latent.obs["cell_type"] = adata.obs[get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
    reference_latent.obs["batch"] = adata.obs[get_from_config(configuration, parameters.CONDITION_KEY)].tolist()
    sc.pp.neighbors(reference_latent, n_neighbors=get_from_config(configuration, parameters.NUMBER_OF_NEIGHBORS))
    sc.tl.leiden(reference_latent)
    sc.tl.umap(reference_latent)
    # no need to show scatterplot during computation
    # sc.pl.umap(reference_latent,
    #           color=['batch', 'cell_type'],
    #           frameon=False,
    #           wspace=0.6,
    #           )

    return reference_latent


def compute_query(pretrained_model, anndata, reference_latent, source_adata, configuration):
    """
    trains the model on a query and saves the result
    :param anndata:
    :return:
    """
    model = sca.models.SCVI.load_query_data(
        anndata,
        'assets/scVI/',
        freeze_dropout=True,
    )
    model.train(
        max_epochs=get_from_config(configuration, parameters.SCVI_QUERY_MAX_EPOCHS),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
        use_gpu=get_from_config(configuration, parameters.USE_GPU)
    )
    tempdir = tempfile.mkdtemp()
    model.save(tempdir, overwrite=True)
    if get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'scvi-model-after-query-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
    utils.delete_file(tempdir + '/model.pt')
    os.removedirs(tempdir)
    if get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            if reference_latent is not None:
                utils.write_latent_csv(reference_latent, key='reference-latent-post-query-training.csv')
            utils.write_adata_to_csv(model, anndata, key='query-adata-post-query-training.csv',
                                     cell_type_key=get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                     condition_key=get_from_config(configuration, parameters.CONDITION_KEY))
            utils.write_adata_to_csv(model, source_adata, key='source-adata-post-query-training.csv',
                                     cell_type_key=get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                     condition_key=get_from_config(configuration, parameters.CONDITION_KEY))
        except Exception as e:
            print(e, file=sys.stderr)

    query_latent = compute_latent(model, anndata, configuration)
    if get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_latent_csv(query_latent, key='query-latent-post-query-training.csv')
            if reference_latent is not None:
                utils.write_combined_csv(reference_latent, query_latent, key='combined-latents-after-query.csv')
        except Exception as e:
            print(e, file=sys.stderr)
    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(query_latent, 'data/figures/query.pdf', color=['batch', 'cell_type'])

    # utils.write_latent_csv(query_latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))
    # utils.write_combined_csv(reference_latent, query_latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))
    # utils.write_adata_to_csv(model, source_adata, key=get_from_config(configuration, parameters.OUTPUT_PATH),
    #                              cell_type_key=get_from_config(configuration, parameters.CELL_TYPE_KEY),
    #                              condition_key=get_from_config(configuration, parameters.CONDITION_KEY))
    utils.write_full_adata_to_csv(model, source_adata, anndata,
                                  key=get_from_config(configuration, parameters.OUTPUT_PATH),
                                  cell_type_key=get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                  condition_key=get_from_config(configuration, parameters.CONDITION_KEY))

    return model


def compute_full_latent(source_adata, target_adata, model, configuration):
    """
    basically just takes to datasets, concatenates them and then computes the latent and saves the result
    :param source_adata:
    :param target_adata:
    :param model:
    :return:
    """
    full_latent = compute_latent(model, source_adata.concatenate(target_adata), configuration)

    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(full_latent, 'data/figures/both.pdf', color=['batch', 'cell_type'])

    both_path = 'both.csv'
    utils.write_latent_csv(full_latent, key='both.csv', filename=both_path)

    return full_latent


def compute_scVI(configuration):
    setup()
    source_adata, target_adata = pre_process_data(configuration)
    model, reference_latent = create_scVI_model(source_adata, target_adata, configuration)
    model = compute_query(model, target_adata, reference_latent, source_adata, configuration)
    # TODO figure out if we need to do this
    # compute_full_latent(source_adata, target_adata, model)
    # model.save('resulting_model', overwrite=True)
    # utils.store_file_in_s3('resulting_model/model.pt', get_from_config(configuration, parameters.RESULTING_MODEL_PATH) + '_new')

