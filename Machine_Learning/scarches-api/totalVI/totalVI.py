import os
import sys
import tempfile
import warnings

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scarches as sca
import torch
from utils import utils, parameters


def setup_modules():
    """
    Set up the warnings filter and the figure parameters
    :return:
    """
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=UserWarning)
    sc.settings.set_figure_params(dpi=200, frameon=False)
    sc.set_figure_params(dpi=200)
    sc.set_figure_params(figsize=(4, 4))
    torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)


def prepare_data(configuration):
    """
    preprocess reference and query dataset
    :param configuration: config
    :return: ref adata, query adata
    """
    # scv.data.pbmcs_10x_cite_seq()
    adata_ref = utils.read_h5ad_file_from_s3(utils.get_from_config(configuration, parameters.REFERENCE_DATA_PATH))
    adata_ref.obs["type"] = "reference"
    # scv.data.dataset_10x("pbmc_10k_v3")
    adata_query = utils.read_h5ad_file_from_s3(utils.get_from_config(configuration, parameters.QUERY_DATA_PATH))
    adata_query.obs["type"] = "query"

    adata_query.obs["batch"] = "PBMC 10k (RNA only)"
    pro_exp = adata_ref.obsm["protein_expression"]  # put matrix of zeros for protein expression (considered missing)
    data = np.zeros((adata_query.n_obs, pro_exp.shape[1]))
    adata_query.obsm["protein_expression"] = pd.DataFrame(columns=pro_exp.columns, index=adata_query.obs_names,
                                                          data=data)

    #concatenate the objects, which intersects the genes properly
    adata_full = anndata.concat([adata_ref, adata_query])

    #split them back up into reference and query
    adata_ref = adata_full[np.logical_or(adata_full.obs.batch == "PBMC5k", adata_full.obs.batch == "PBMC10k")].copy()
    adata_query = adata_full[adata_full.obs.batch == "PBMC 10k (RNA only)"].copy()

    # run gene selection on the reference
    sc.pp.highly_variable_genes(
        adata_ref,
        n_top_genes=4000,
        flavor="seurat_v3",
        batch_key="batch",
        subset=True,
    )
    #use selected genes for the query dataset
    adata_query = adata_query[:, adata_ref.var_names].copy()
    return adata_ref, adata_query


def train_model(adata_ref, configuration):
    """
    Create TOTALVI model and train it on reference dataset
    :param adata_ref: reference dataset
    :param configuration: config
    :return: trained totalVI model
    """
    sca.models.TOTALVI.setup_anndata(
        adata_ref,
        batch_key="batch",
        protein_expression_obsm_key="protein_expression"
    )
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none",
    )
    if utils.get_from_config(configuration, parameters.USE_PRETRAINED_TOTALVI_MODEL):
        vae_ref = sca.models.TOTALVI.load(adata=adata_ref, dir_path='assets/totalVI/' + str(
            utils.get_from_config(configuration, parameters.ATLAS)) + '/')
    else:
        vae_ref = sca.models.TOTALVI(adata_ref, **arches_params)
        vae_ref.train(utils.get_from_config(configuration, parameters.TOTALVI_MAX_EPOCHS_1),
                      use_gpu=utils.get_from_config(configuration, parameters.USE_GPU))
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.write_adata_to_csv(vae_ref, adata_ref, key='totalvi-source-adata-post-first-training.csv')
            except Exception as e:
                print(e, file=sys.stderr)
        tempdir = tempfile.mkdtemp()
        vae_ref.save(tempdir, overwrite=True)
        if utils.get_from_config(configuration, parameters.DEV_DEBUG):
            try:
                utils.store_file_in_s3(tempdir + '/model.pt', 'totalvi-model-after-first-training.pt')
            except Exception as e:
                print(e, file=sys.stderr)
        utils.delete_file(tempdir + '/model.pt')
        os.removedirs(tempdir)
    return vae_ref


def visualize_and_store_as_pdf(filename, adata_ref, **config):
    """
    visualize reference dataset and save visualization with new file name
    :param filename: our file name
    :param adata_ref: reference data
    :param config: config
    :return:
    """
    sc.pl.umap(adata_ref, **config)
    os.rename('figures/umap.pdf', 'figures/' + filename)
    # os.rename(args.pdfpath + 'umap.pdf', args.pdfpath + filename)


def visualize_RNA_data(model, adata_ref, configuration):
    """
    Save Latent representation and visualize RNA data
    :param model: totalVI model
    :param adata_ref: reference dataset
    :param configuration: config
    :return:
    """
    adata_ref.obsm["X_totalVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata_ref, use_rep="X_totalVI")
    sc.tl.umap(adata_ref, min_dist=0.4)
    # utils.write_latent_csv(adata_ref, key=utils.get_from_config(configuration, parameters.OUTPUT_PATH))
    # utils.write_full_adata_to_csv(model, adata_ref, sc.AnnData(model.get_latent_representation()), key=utils.get_from_config(configuration, parameters.OUTPUT_PATH),
    #                              cell_type_key=utils.get_from_config(configuration, parameters.CELL_TYPE_KEY),
    #                              condition_key=utils.get_from_config(configuration, parameters.CONDITION_KEY))
    if utils.get_from_config(configuration, parameters.DEBUG):
        visualize_and_store_as_pdf("firstumap.pdf",
                                   adata_ref,
                                   color=["batch"],
                                   frameon=False,
                                   ncols=1,
                                   title="Reference",
                                   save=True)


def surgery(adata_query, configuration):
    """
    Perform surgery on reference model and train on query dataset without protein data
    :param adata_query: query dataset
    :param configuration: config
    :return: trained model on query dataset
    """
    dir_path = 'assets/totalVI/' + utils.get_from_config(configuration, parameters.ATLAS) + '/'
    vae_q = sca.models.TOTALVI.load_query_data(
        adata_query,
        dir_path,
        freeze_expression=True
    )
    vae_q.train(int(utils.get_from_config(configuration, parameters.TOTALVI_MAX_EPOCHS_2)),
                plan_kwargs=dict(weight_decay=0.0), use_gpu=utils.get_from_config(configuration, parameters.USE_GPU))
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.write_adata_to_csv(vae_q, adata_query, key='query-adata-post-second-training.csv')
        except Exception as e:
            print(e, file=sys.stderr)
    tempdir = tempfile.mkdtemp()
    vae_q.save(tempdir, overwrite=True)
    if utils.get_from_config(configuration, parameters.DEV_DEBUG):
        try:
            utils.store_file_in_s3(tempdir + '/model.pt', 'totalvi-model-after-query-training.pt')
        except Exception as e:
            print(e, file=sys.stderr)
    utils.delete_file(tempdir + '/model.pt')
    os.removedirs(tempdir)
    adata_query.obsm["X_totalVI"] = vae_q.get_latent_representation()
    sc.pp.neighbors(adata_query, use_rep="X_totalVI")
    sc.tl.umap(adata_query, min_dist=0.4)
    return vae_q


def impute_proteins(vae_q, adata_query, configuration):
    """
    Impute protein data for the query dataset and visualize
    :param vae_q: trained model
    :param adata_query: query data
    :param configuration: config
    :return:
    """
    #Impute the proteins that were observed in the reference,
    # using the transform_batch parameter
    _, imputed_proteins = vae_q.get_normalized_expression(
        adata_query,
        n_samples=25,
        return_mean=True,
        transform_batch=["PBMC10k", "PBMC5k"],
    )
    adata_query.obs = pd.concat([adata_query.obs, imputed_proteins], axis=1)
    if utils.get_from_config(configuration, parameters.DEBUG):
        visualize_and_store_as_pdf("secondumap.pdf",
                                   adata_query,
                                   color=imputed_proteins.columns,
                                   frameon=False,
                                   ncols=3,
                                   save=True)


def latent_ref_representation(adata_query, adata_ref, vae_q):
    """
    Get latent representation of reference + query dataset
    and compute UMAP
    :param adata_query: query data
    :param adata_ref: reference data
    :param vae_q: trained model
    :return: concatenated query data, imputed proteins
    """
    adata_full_new = adata_query.concatenate(adata_ref, batch_key="none")
    adata_full_new.obsm["X_totalVI"] = vae_q.get_latent_representation(adata=adata_full_new)
    sc.pp.neighbors(adata_full_new, use_rep="X_totalVI")
    sc.tl.umap(adata_full_new, min_dist=0.3)
    _, imputed_proteins_all = vae_q.get_normalized_expression(
        adata_full_new,
        n_samples=25,
        return_mean=True,
        transform_batch=["PBMC10k", "PBMC5k"],
    )
    for i, p in enumerate(imputed_proteins_all.columns):
        adata_full_new.obs[p] = imputed_proteins_all[p].to_numpy().copy()
    return adata_full_new, imputed_proteins_all


def compute_final_umaps(adata_full_new, imputed_proteins_all, configuration):
    """
     save results and write latent into csv file
    :param adata_full_new: concatenated query data
    :param imputed_proteins_all:
    :param configuration:
    :return:
    """
    perm_inds = np.random.permutation(np.arange(adata_full_new.n_obs))
    utils.write_latent_csv(adata_full_new[perm_inds], key=utils.get_from_config(configuration, parameters.OUTPUT_PATH))
    if utils.get_from_config(configuration, parameters.DEBUG):
        visualize_and_store_as_pdf("thirdumap.pdf",
                                   adata_full_new[perm_inds],
                                   color=["batch"],
                                   frameon=False,
                                   ncols=1,
                                   title="Reference and query",
                                   save=True
                                   )
        """
        ax = visualize_and_store_as_pdf("forthumap.pdf",
                                        adata_full_new,
                                        color="batch",
                                        groups=["PBMC 10k (RNA only)"],
                                        frameon=False,
                                        ncols=1,
                                        title="Reference and query",
                                        alpha=0.4,
                                        save=True
                                        )
        visualize_and_store_as_pdf("fifthumap.pdf",
                                   adata_full_new,
                                   color=imputed_proteins_all.columns,
                                   frameon=False,
                                   ncols=3,
                                   vmax="p99",
                                   save=True
                                   )
        """


def computeTotalVI(configuration):
    # this is not needed when running inside api
    # global args
    # args = parse_args()
    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.DEBUG if args.debug else logging.INFO)
    # if not args.example and not (args.ref or args.query)(
    #        exists(utils.get_from_config(configuration, parameters.REFERENCE_DATA_PATH)) or exists(
    #                utils.get_from_config(configuration, parameters.QUERY_DATA_PATH))):  
    #    logger.error("file path to 'ref' and 'query' can't be empty if the argument 'example' is set to false")
    #    exit()
    """
    process reference and query dataset with totalVI model
    :param configuration: config
    :return:
    """
    setup_modules()
    adata_ref, adata_query = prepare_data(configuration)
    model_ref = train_model(adata_ref, configuration)
    visualize_RNA_data(model_ref, adata_ref, configuration)
    vae_q = surgery(adata_query, configuration)
    impute_proteins(vae_q, adata_query, configuration)
    adata_full_new, imputed_proteins_all = latent_ref_representation(adata_query, adata_ref, vae_q)
    compute_final_umaps(adata_full_new, imputed_proteins_all, configuration)
