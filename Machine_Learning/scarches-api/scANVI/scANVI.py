import warnings
from scVI.scVI import create_scVI_model, set_config

import scanpy
import scarches
from scarches.dataset.trvae.data_handling import remove_sparsity
from matplotlib import pyplot as plt
import numpy as np
import torch
import utils
import logging

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


def pre_process_data():
    source_adata = scanpy.read(get_from_config('reference_dataset'))
    target_adata = scanpy.read(get_from_config('query_dataset'))
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

    # model.save(get_from_config('ref_path'), overwrite=True)     # würde das speichern wo anders machen, muss ja an verschiedenen Orten gespeichert werden

    return reference_latent


def predict(model, latent):
    latent.obs['predictions'] = model.predict()
    print("Acc: {}".format(np.mean(latent.obs.predictions == latent.obs.cell_type)))
    return latent


def surgery(anndata):
    model = scarches.models.SCANVI.load_query_data(
        anndata,
        get_from_config('model_path'),  # ist das der richtige Pfad? Ist doch dann schon einmal trainiert?
        freeze_dropout=True,
    )

    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []

    print("Labelled Indices: ", len(model._labeled_indices))
    print("Unlabelled Indices: ", len(model._unlabeled_indices))

    model.train(
        max_epochs=get_from_config('scanvi_max_epochs'),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )

    surgery_latent = get_latent(model, anndata)

    if get_from_config('debug'):
        utils.save_umap_as_pdf(surgery_latent, 'figures/surgery.pdf', color=['batch', 'cell_type'])
        utils.write_latent_csv(surgery_latent, filename='/tmp/surgery.csv')

    model.save(get_from_config('surgery_path'), overwrite=True)

    return model, surgery_latent


def query(anndata):
    model = scarches.models.SCANVI.load_query_data(
        anndata,
        get_from_config('model_path'),
        # ist das der richtige Pfad? Das wäre ja dann das gerade von scVI convertierte, hier gabs nen Fehler von wegen falsche Klasse (also weil das convertierte scanVI wahrscheinlich nicht gespeichert wurde)
        freeze_dropout=True,
        # habs jetzt mal unten gespeichert bevor query aufgerufen wird, schau es dir aber nochmal an
    )

    model._unlabeled_indices = np.arange(anndata.n_obs)
    model._labeled_indices = []

    if get_from_config('debug'):
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", len(model._unlabeled_indices))

    model.train(
        max_epochs=get_from_config('scanvi_max_epochs'),
        plan_kwargs=dict(weight_decay=0.0),
        check_val_every_n_epoch=10,
    )

    query_latent = get_latent(model, anndata)

    if get_from_config('debug'):
        utils.save_umap_as_pdf(query_latent, 'figures/query.pdf', color=['batch', 'cell_type'])
        utils.write_latent_csv(query_latent, filename='/tmp/query.csv')

        # Muss man das Model dann nicht abspeichern? Oder ist das dann nicht mehr das pre-trained?

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

    if get_from_config('debug'):
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

    if get_from_config('debug'):
        utils.save_umap_as_pdf(latent, 'figures/compare.pdf', color=["predictions", "cell_type"])
        utils.write_latent_csv(latent, filename='/tmp/compare.csv')


def compute_scANVI(configP):
    global config
    config = configP

    if get_from_config('debug'):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

    setup_modules()

    source_adata, target_adata = pre_process_data()
    #                                                   if not get_from_config('pre_trained_scANVI'): bin mir bei dir nicht sicher, was alles zu was gehöhrt
    set_config(configP)  # sets the config in scVI

    if get_from_config('debug'):
        logger.debug(source_adata)
        logger.debug(target_adata)

    vae = create_scVI_model(source_adata,
                            target_adata)  # kann man sehen ob es schon ein scVI oder ein scANVI model gibt?

    if get_from_config('debug'):
        logger.debug(source_adata)
        logger.debug(target_adata)

    # if args.train:
    # vae.train(max_epochs=get_from_config('scvi_max_epochs'))

    # setup_anndata_for_scanvi(source_adata)

    scanvi = get_scanvi_from_scvi_model(vae)

    if get_from_config('debug'):
        print("Labelled Indices: ", len(scanvi._labeled_indices))
        print("Unlabelled Indices: ", len(scanvi._unlabeled_indices))

    scanvi.train(max_epochs=get_from_config('scanvi_max_epochs'))

    reference_latent = get_latent(scanvi, source_adata)

    if get_from_config('debug'):
        utils.save_umap_as_pdf(reference_latent, 'figures/reference.pdf', color=['batch', 'cell_type'])

    reference_latent = predict(scanvi, reference_latent)

    scanvi.save(get_from_config('ref_path'),
                overwrite=True)  # ich bin mir nicht sicher, wann dein Model in welchem Stadium ist

    model_query, query_latent = query(target_adata)
    model = model_query

    if predict:
        predict_latent(predict(model_query, query_latent))

    model_surgery, surgery_latent = surgery(target_adata)

    if predict:
        predict_latent(predict(model_surgery, surgery_latent))

    full_latent = None

    if get_from_config('both'):
        full_latent = both_adata(source_adata, target_adata)

    if get_from_config('compare'):
        if full_latent is None:
            full_latent = both_adata(source_adata, target_adata)
        if model is None:
            model_query, query_latent = query(target_adata)
            model = model_query
        compare_adata(model, source_adata, target_adata, full_latent)
