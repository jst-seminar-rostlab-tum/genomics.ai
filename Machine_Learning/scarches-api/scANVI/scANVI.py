import os
import warnings
from scVI import scVI
import scanpy
import scarches
from scarches.dataset.trvae.data_handling import remove_sparsity
from matplotlib import pyplot as plt
import numpy as np
import torch
from utils import utils, parameters
import logging
import tempfile
import sys


def get_from_config(configuration, key):
    if key in configuration:
        return configuration[key]
    return None


def setup_modules():
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)

    scanpy.settings.set_figure_params(dpi=200, frameon=False)
    scanpy.set_figure_params(dpi=200)
    scanpy.set_figure_params(figsize=(4, 4))

    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)


def pre_process_data(configuration):
    source_adata = utils.read_h5ad_file_from_s3(get_from_config(configuration, parameters.REFERENCE_DATA_PATH))
    target_adata = utils.read_h5ad_file_from_s3(get_from_config(configuration, parameters.QUERY_DATA_PATH))
    source_adata = remove_sparsity(source_adata)
    target_adata = remove_sparsity(target_adata)

    return source_adata, target_adata


def setup_anndata_for_scanvi(anndata, configuration):
    scarches.models.SCANVI.setup_anndata(anndata,
                                         unlabeled_category=get_from_config(configuration, parameters.UNLABELED_KEY),
                                         batch_key=get_from_config(configuration, parameters.CONDITION_KEY),
                                         labels_key=get_from_config(configuration, parameters.CELL_TYPE_KEY))


def get_scanvi_from_scvi_model(scvi_model, configuration):
    return scarches.models.SCANVI.from_scvi_model(scvi_model, get_from_config(configuration, parameters.UNLABELED_KEY))


def get_latent(model, adata, configuration):
    reference_latent = scanpy.AnnData(model.get_latent_representation())
    reference_latent.obs["cell_type"] = adata.obs[get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
    reference_latent.obs["batch"] = adata.obs[get_from_config(configuration, parameters.CONDITION_KEY)].tolist()

    scanpy.pp.neighbors(reference_latent, n_neighbors=get_from_config(configuration, parameters.NUMBER_OF_NEIGHBORS))
    scanpy.tl.leiden(reference_latent)
    scanpy.tl.umap(reference_latent)

    return reference_latent


def predict(model, latent):
    latent.obs['predictions'] = model.predict()
    print("Acc: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
    return latent


def surgery(reference_latent, source_adata, anndata, configuration):
    model = scarches.models.SCANVI.load_query_data(
        anndata,
        get_from_config(configuration, parameters.PRETRAINED_MODEL_PATH),
        # ist das der richtige Pfad? Ist doch dann schon einmal trainiert?
        freeze_dropout=True,
    )

    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []

    print("Labelled Indices: ", len(model._labeled_indices))
    print("Unlabelled Indices: ", len(model._unlabeled_indices))

    model.train(
        max_epochs=get_from_config(configuration, parameters.SCANVI_MAX_EPOCHS),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
        use_gpu=get_from_config(configuration, parameters.USE_GPU)
    )

    surgery_latent = get_latent(model, anndata, configuration)

    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(surgery_latent, 'figures/surgery.pdf', color=['batch', 'cell_type'])

    # utils.write_combined_csv(reference_latent, surgery_latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))
    utils.write_full_adata_to_csv(model, source_adata, anndata,
                                  key=get_from_config(configuration, parameters.OUTPUT_PATH),
                                  cell_type_key=get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                  condition_key=get_from_config(configuration, parameters.CONDITION_KEY))

    model.save('scvi_model', overwrite=True)  # TODO check if we need this, for now, we delete it
    utils.delete_file('scvi_model/model.pt')
    os.rmdir('scvi_model')

    return model, surgery_latent


def query(pretrained_model, reference_latent, anndata, source_adata, configuration):
    model = pretrained_model

    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []

    if get_from_config(configuration, parameters.DEBUG):
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", len(model._unlabeled_indices))

    model.train(
        max_epochs=get_from_config(configuration, parameters.SCANVI_MAX_EPOCHS_QUERY),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
        use_gpu=get_from_config(configuration, parameters.USE_GPU)
    )
    tempdir = tempfile.mkdtemp()
    model.save(tempdir, overwrite=True)
    if get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_adata_to_csv(model, 'scanvi-query-latent-after-query-training.csv')
        except Exception as e:
            print(e, file=sys.stderr)
    if get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'scanvi-model-after-query-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
    utils.delete_file(tempdir + '/model.pt')
    os.removedirs(tempdir)
    query_latent = get_latent(model, anndata, configuration)

    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(query_latent, 'figures/query.pdf', color=['batch', 'cell_type'])

    # utils.write_combined_csv(reference_latent, query_latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))
    utils.write_full_adata_to_csv(model, source_adata, anndata,
                                  key=get_from_config(configuration, parameters.OUTPUT_PATH),
                                  cell_type_key=get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                  condition_key=get_from_config(configuration, parameters.CONDITION_KEY))

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


def both_adata(source_adata, target_adata, configuration):
    adata_full = source_adata.concatenate(target_adata)
    full_latent = scanpy.AnnData(scarches.models.SCANVI.get_latent_representation(adata=adata_full))
    full_latent.obs['cell_type'] = adata_full.obs[get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
    full_latent.obs['batch'] = adata_full.obs[get_from_config(configuration, parameters.CONDITION_KEY)].tolist()

    scanpy.pp.neighbors(full_latent)
    scanpy.tl.leiden(full_latent)
    scanpy.tl.umap(full_latent)

    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(full_latent, 'figures/both.pdf', color=['batch', 'cell_type'])

    utils.write_latent_csv(full_latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))

    return full_latent


def compare_adata(model, source_adata, target_adata, latent, configuration):
    adata_full = source_adata.concatenate(target_adata)
    latent.obs['predictions'] = model.predict(adata=adata_full)

    print("Acc: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
    scanpy.pp.neighbors(latent)
    scanpy.tl.leiden(latent)
    scanpy.tl.umap(latent)

    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(latent, 'figures/compare.pdf', color=["predictions", "cell_type"])

    utils.write_latent_csv(latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))


def create_model(source_adata, target_adata, configuration):
    if get_from_config(configuration, parameters.USE_PRETRAINED_SCANVI_MODEL):
        return scarches.models.SCANVI.load_query_data(
            target_adata,
            'assets/scANVI/',
            freeze_dropout=True,
        ), None

    scvi_model, _ = scVI.create_scVI_model(source_adata, target_adata, configuration)
    scanvi = get_scanvi_from_scvi_model(scvi_model, configuration)

    if get_from_config(configuration, parameters.DEBUG):
        print("Labelled Indices: ", len(scanvi._labeled_indices))
        print("Unlabelled Indices: ", len(scanvi._unlabeled_indices))

    scanvi.train(max_epochs=get_from_config(configuration, parameters.SCANVI_MAX_EPOCHS),
                 use_gpu=get_from_config(configuration, parameters.USE_GPU))
    tempdir = tempfile.mkdtemp()
    scanvi.save(tempdir, overwrite=True)
    if get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_adata_to_csv(scanvi, 'scanvi-reference-latent-after-from-scvi-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'scanvi-model-after-first-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)

    utils.delete_file(tempdir + '/model.pt')
    os.removedirs(tempdir)

    reference_latent = get_latent(scanvi, source_adata, configuration)

    if get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(reference_latent, 'figures/reference.pdf', color=['batch', 'cell_type'])

    # TODO check if needed
    # reference_latent = predict(scanvi, reference_latent)
    return scanvi, reference_latent


def compute_scANVI(configuration):
    if get_from_config(configuration, parameters.DEBUG):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

    setup_modules()

    source_adata, target_adata = pre_process_data(configuration)

    scanvi, reference_latent = create_model(source_adata, target_adata, configuration)

    model_query, query_latent = query(scanvi, reference_latent, target_adata, source_adata, configuration)
    model = model_query

    if get_from_config(configuration, parameters.SCANVI_PREDICT_CELLTYPES):
        predict_latent(predict(model_query, query_latent))

    #if get_from_config(configuration, parameters.SCANVI_DO_SURGERY):
    #    model_surgery, surgery_latent = surgery(reference_latent, target_adata, configuration)
    #
    #    if get_from_config(configuration, parameters.SCANVI_PREDICT_CELLTYPES):
    #        predict_latent(predict(model_surgery, surgery_latent))

    full_latent = None

    # TODO check if needed at all
    #if get_from_config(configuration, parameters.SCANVI_COMPARE_REFERENCE_AND_QUERY):
    #    full_latent = both_adata(source_adata, target_adata, configuration)

    #if get_from_config(configuration, parameters.SCANVI_COMPARE_OBSERVED_AND_PREDICTED_CELLTYPES):
    #    if full_latent is None:
    #        full_latent = both_adata(source_adata, target_adata, configuration)
    #    if model is None:
    #        model_query, query_latent = query(None, reference_latent, target_adata, configuration)
    #        model = model_query
    #    compare_adata(model, source_adata, target_adata, full_latent, configuration)
