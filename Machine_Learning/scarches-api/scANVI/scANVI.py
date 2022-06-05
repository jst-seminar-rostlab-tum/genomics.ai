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
import scvi


def setup_modules():
    """
    Set up the warnings filter and the figure parameters
    :return:
    """
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)

    scanpy.settings.set_figure_params(dpi=200, frameon=False)
    scanpy.set_figure_params(dpi=200)
    scanpy.set_figure_params(figsize=(4, 4))

    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)


def setup_anndata_for_scanvi(anndata, configuration):
    """
    Preprocess reference dataset
    :param anndata: reference dataset /source adata
    :param configuration: config
    :return:
    """
    scarches.models.SCANVI.setup_anndata(anndata,
                                         batch_key='dataset',
                                         labels_key='scanvi_label',
                                         unlabeled_category='unlabeled')


def get_scanvi_from_scvi_model(scvi_model, configuration):
    """
    Create the scANVI model instance
    :param scvi_model: the scVI model
    :param configuration: config
    :return: scANVI model
    """
    return scarches.models.SCANVI.from_scvi_model(scvi_model,
                                                  utils.get_from_config(configuration, parameters.UNLABELED_KEY))


def get_latent(model, adata, configuration):
    """
    Create anndata file of latent representation and compute UMAP
    :param model: the created scANVI model
    :param adata: reference dataset / source adata
    :param configuration: config
    :return: latent representation
    """
    # add obs to reference_latent
    reference_latent = scanpy.AnnData(model.get_latent_representation())
    reference_latent.obs["cell_type"] = adata.obs[
        utils.get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
    reference_latent.obs["batch"] = adata.obs[utils.get_from_config(configuration, parameters.CONDITION_KEY)].tolist()

    scanpy.pp.neighbors(reference_latent,
                        n_neighbors=utils.get_from_config(configuration, parameters.NUMBER_OF_NEIGHBORS))
    scanpy.tl.leiden(reference_latent)
    scanpy.tl.umap(reference_latent)

    return reference_latent


def predict(model, latent):
    """
    predict on the latent and compute the accuracy of the predicted value
    :param model: scANVI model
    :param latent: reference latent
    :return: latent with obs "predicted"
    """
    latent.obs['predicted'] = model.predict()
    print("Acc: {}".format(np.mean(latent.obs.predicted == latent.obs.cell_type)))
    return latent


def surgery(reference_latent, source_adata, anndata, configuration):
    """
    Perform surgery on reference model
    :param reference_latent: reference latent
    :param source_adata: reference dataset
    :param anndata: query dataset
    :param configuration: config
    :return: trained model, surgery latent
    """
    model = scarches.models.SCANVI.load_query_data(
        anndata,
        utils.get_from_config(configuration, parameters.PRETRAINED_MODEL_PATH),
        # ist das der richtige Pfad? Ist doch dann schon einmal trainiert?
        freeze_dropout=True,
    )

    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []

    print("Labelled Indices: ", len(model._labeled_indices))
    print("Unlabelled Indices: ", len(model._unlabeled_indices))

    model.train(
        max_epochs=utils.get_from_config(configuration, parameters.SCANVI_MAX_EPOCHS),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
        use_gpu=utils.get_from_config(configuration, parameters.USE_GPU)
    )

    surgery_latent = get_latent(model, anndata, configuration)

    if utils.get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(surgery_latent, 'figures/surgery.pdf', color=['batch', 'cell_type'])

    # utils.write_combined_csv(reference_latent, surgery_latent, key=utils.get_from_config(configuration, parameters.OUTPUT_PATH))
    utils.write_full_adata_to_csv(model, source_adata, anndata,
                                  key=utils.get_from_config(configuration, parameters.OUTPUT_PATH),
                                  cell_type_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                  condition_key=utils.get_from_config(configuration, parameters.CONDITION_KEY),
                                  predictScanvi=True, configuration=configuration)

    model.save('scvi_model', overwrite=True)
    utils.delete_file('scvi_model/model.pt')
    os.rmdir('scvi_model')

    return model, surgery_latent


def query(pretrained_model, reference_latent, anndata, source_adata, configuration):
    """
    Perform surgery on reference model and train on query dataset
    :param pretrained_model: pretrained model
    :param reference_latent: reference latent
    :param anndata: target adata / query dataset
    :param source_adata: reference dataset
    :param configuration: config
    :return: trained model, query latent
    """
    model = scarches.models.SCANVI.load_query_data(
        anndata,
        'assets/scANVI/' + str(utils.get_from_config(configuration, parameters.ATLAS)) + '/',
        freeze_dropout=True,
    )

    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []

    if utils.get_from_config(configuration, parameters.DEBUG):
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", len(model._unlabeled_indices))

#TODO: HARDCODING for human lung cell atlas -------------------------------------
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
            max_epochs=utils.get_from_config(configuration, parameters.SCANVI_MAX_EPOCHS_QUERY),
            plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10,
            use_gpu=utils.get_from_config(configuration, parameters.USE_GPU)
        )
    tempdir = tempfile.mkdtemp()
    model.save(tempdir, overwrite=True)
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_adata_to_csv(model, 'scanvi-query-latent-after-query-training.csv')
        except Exception as e:
            print(e, file=sys.stderr)
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'scanvi-model-after-query-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
    utils.delete_file(tempdir + '/model.pt')
    os.removedirs(tempdir)

    # add obs to query_latent
    query_latent = get_latent(model, anndata, configuration)
    query_latent.obs['cell_type'] = anndata.obs[utils.get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
    query_latent.obs['batch'] = anndata.obs[utils.get_from_config(configuration, parameters.CONDITION_KEY)].tolist()
    scanpy.pp.neighbors(query_latent, n_neighbors=utils.get_from_config(configuration, parameters.NUMBER_OF_NEIGHBORS))
    scanpy.tl.leiden(query_latent)
    scanpy.tl.umap(query_latent)

    if utils.get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(query_latent, 'figures/query.pdf', color=['batch', 'cell_type'])
    utils.write_full_adata_to_csv(model, source_adata, anndata,
                                  key=utils.get_from_config(configuration, parameters.OUTPUT_PATH),
                                  cell_type_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY),
                                  condition_key=utils.get_from_config(configuration, parameters.CONDITION_KEY),
                                  predictScanvi=True, configuration=configuration)
    return model, query_latent


def predict_latent(model, latent):
    """
    Compute Accuracy of model for query dataset
    compare predicted and observed cell types
    :param model: scANVI model
    :param latent: query latent
    :return:
    """
    latent.obs['predicted'] = model.predict()
    print("Acc: {}".format(np.mean(latent.obs.predicted == latent.obs.cell_type)))

    df = latent.obs.groupby(["cell_type", "predicted"]).size().unstack(fill_value=0)
    norm_df = df / df.sum(axis=0)

    figure = plt.figure(figsize=(8, 8), frameon=False)
    _ = plt.grid(False)
    _ = plt.pcolor(norm_df)
    _ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
    _ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
    plt.xlabel("Predicted")
    plt.ylabel("Observed")
    figure.savefig('predict.png')


# def both_adata(source_adata, target_adata, configuration):
#     adata_full = source_adata.concatenate(target_adata)
#     full_latent = scanpy.AnnData(scarches.models.SCANVI.get_latent_representation(adata=adata_full))
#     full_latent.obs['cell_type'] = adata_full.obs[
#         utils.get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
#     full_latent.obs['batch'] = adata_full.obs[utils.get_from_config(configuration, parameters.CONDITION_KEY)].tolist()

#     scanpy.pp.neighbors(full_latent)
#     scanpy.tl.leiden(full_latent)
#     scanpy.tl.umap(full_latent)

#     full_latent.obs['predicted'] = 'predicted'

#     if utils.get_from_config(configuration, parameters.DEBUG):
#         utils.save_umap_as_pdf(full_latent, 'figures/both.pdf', color=['batch', 'cell_type'])

#     utils.write_latent_csv(full_latent, key=utils.get_from_config(configuration, parameters.OUTPUT_PATH))

#     return full_latent


# def compare_adata(model, source_adata, target_adata, configuration):
#     adata_full = source_adata.concatenate(target_adata)
#     full_latent = scanpy.AnnData(scarches.models.SCANVI.get_latent_representation(adata=adata_full))
#     full_latent.obs['cell_type'] = adata_full.obs[
#         utils.get_from_config(configuration, parameters.CELL_TYPE_KEY)].tolist()
#     full_latent.obs['batch'] = adata_full.obs[utils.get_from_config(configuration, parameters.CONDITION_KEY)].tolist()

#     scanpy.pp.neighbors(full_latent)
#     scanpy.tl.leiden(full_latent)
#     scanpy.tl.umap(full_latent)
#
#     full_latent.obs['predictions'] = 'predicted'
#
#
#     latent.obs['predictions'] = model.predict(adata=adata_full)
#     print("Acc_compare: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
#     scanpy.pp.neighbors(latent)
#     scanpy.tl.leiden(latent)
#     scanpy.tl.umap(latent)
#
#     if get_from_config(configuration, parameters.DEBUG):
#         utils.save_umap_as_pdf(latent, 'figures/compare.pdf', color=["predictions", "cell_type"])
#
#     utils.write_latent_csv(latent, key=get_from_config(configuration, parameters.OUTPUT_PATH))

#     full_latent.obs['predicted'] = model.predict(adata=adata_full)
#     print("Acc_compare: {}".format(np.mean(full_latent.obs.predicted == full_latent.obs.cell_type)))
#     scanpy.pp.neighbors(full_latent)
#     scanpy.tl.leiden(full_latent)
#     scanpy.tl.umap(full_latent)

#     if utils.get_from_config(configuration, parameters.DEBUG):
#         utils.save_umap_as_pdf(full_latent, 'figures/compare.pdf', color=["predicted", "cell_type"])

#     utils.write_latent_csv(full_latent, key=utils.get_from_config(configuration, parameters.OUTPUT_PATH))


def create_model(source_adata, target_adata, configuration):
    """
    - compute scANVI model and train it on reference dataset
    - compute the accuracy of the learned classifier
    - save the result and write into csv file to the s3

    :param source_adata: reference dataset
    :param target_adata: query dataset
    :param configuration: config
    :return: scANVI model, reference latent
    """
    if utils.get_from_config(configuration, parameters.USE_PRETRAINED_SCANVI_MODEL):
        path = 'assets/scANVI/' + str(utils.get_from_config(configuration, parameters.ATLAS)) + '/'
        return scarches.models.SCANVI.load_query_data(
            target_adata,
            path,
            freeze_dropout=True,
        ), None

    scvi_model, _ = scVI.create_scVI_model(source_adata, target_adata, configuration)
    scanvi = get_scanvi_from_scvi_model(scvi_model, configuration)

    if utils.get_from_config(configuration, parameters.DEBUG):
        print("Labelled Indices: ", len(scanvi._labeled_indices))
        print("Unlabelled Indices: ", len(scanvi._unlabeled_indices))

    scanvi.train(max_epochs=utils.get_from_config(configuration, parameters.SCANVI_MAX_EPOCHS),
                 use_gpu=utils.get_from_config(configuration, parameters.USE_GPU))
    tempdir = tempfile.mkdtemp()
    scanvi.save(tempdir, overwrite=True, save_anndata=True)
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_adata_to_csv(scanvi, 'scanvi-reference-latent-after-from-scvi-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'scanvi-model-after-first-training.pt')
            utils.store_file_in_s3(tempdir + '/adata.h5ad', 'scanvi-adata-after-first-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)

    utils.delete_file(tempdir + '/model.pt')
    utils.delete_file(tempdir + '/adata.h5ad')
    os.removedirs(tempdir)

    reference_latent = get_latent(scanvi, source_adata, configuration)

    if utils.get_from_config(configuration, parameters.DEBUG):
        utils.save_umap_as_pdf(reference_latent, 'figures/reference.pdf', color=['batch', 'cell_type'])

    reference_latent = predict(scanvi, reference_latent)
    return scanvi, reference_latent


def compute_scANVI(configuration):
    """
    process reference and query dataset with scANVI model
    :param configuration: config
    :return:
    """
    if utils.get_from_config(configuration, parameters.DEBUG):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

    setup_modules()

    source_adata, target_adata = utils.pre_process_data(configuration)

    scanvi, reference_latent = create_model(source_adata, target_adata, configuration)

    model_query, query_latent = query(scanvi, reference_latent, target_adata, source_adata, configuration)
