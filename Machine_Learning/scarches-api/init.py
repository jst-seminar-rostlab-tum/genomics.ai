import sys

from scANVI.scANVI import compute_scANVI
from scVI.scVI import compute_scVI
from totalVI.totalVI import computeTotalVI
from utils import utils, parameters


def default_config():
    """
    returns the default config combined for all the models
    :return: dict containing all the default values
    """
    return {
        parameters.MODEL: 'scANVI',

        parameters.REFERENCE_DATA_PATH: 'data/ref/source_data.h5ad',
        parameters.QUERY_DATA_PATH: 'data/query/target_data.h5ad',
        parameters.RESULTING_MODEL_PATH: 'data/model',
        parameters.OUTPUT_PATH: 'query.tsv',

        parameters.USE_PRETRAINED_SCVI_MODEL: False,
        parameters.USE_PRETRAINED_TOTALVI_MODEL: False,
        parameters.USE_PRETRAINED_SCANVI_MODEL: False,

        # scANVI stuff
        parameters.SCANVI_COMPARE_REFERENCE_AND_QUERY: False,
        parameters.SCANVI_COMPARE_OBSERVED_AND_PREDICTED_CELLTYPES: False,
        parameters.SCANVI_PREDICT_CELLTYPES: False,

        parameters.CONDITION_KEY: 'study',
        parameters.CELL_TYPE_KEY: 'cell_type',
        parameters.TARGET_CONDITIONS: ['Pancreas CelSeq2', 'Pancreas SS2'],
        parameters.PRETRAINED_MODEL_PATH: './ref_model/',
        parameters.NUMBER_OF_LAYERS: 2,
        parameters.ENCODE_COVARIATES: True,
        parameters.DEEPLY_INJECT_COVARIATES: False,
        parameters.USE_LAYER_NORM: 'both',
        parameters.USE_BATCH_NORM: 'none',
        parameters.UNLABELED_KEY: 'Unknown',
        parameters.SCANVI_MAX_EPOCHS: 20,
        parameters.SCVI_MAX_EPOCHS: 20,
        parameters.NUMBER_OF_NEIGHBORS: 8,
        parameters.MAX_EPOCHS: 100,
        parameters.UNWANTED_LABELS: ['leiden'],
        parameters.DEBUG: False,
        parameters.RUN_ASYNCHRONOUSLY: False,
        parameters.ATTRIBUTES: None,

        # totalVI stuff
        parameters.TOTALVI_MAX_EPOCHS_1: 1,  # 400
        parameters.TOTALVI_MAX_EPOCHS_2: 1,  # 200
        parameters.SCANVI_DO_SURGERY: False,
    }


def get_from_config(configuration, key):
    if key in configuration:
        return configuration[key]
    return None


def merge_configs(user_config):
    """
    overwrites the default config with the input from the rest api
    :param user_config:
    :return: dict
    """
    return default_config() | user_config


# def query(reference_dataset, query_dataset, model_path, surgery_path,  model_type):
def query(user_config):
    configuration = merge_configs(user_config)
    model = get_from_config(configuration, parameters.MODEL)
    attributes = None
    if model == 'scVI':
        attributes = compute_scVI(configuration)
    elif model == 'scANVI':
        attributes = compute_scANVI(configuration)
    elif model == 'totalVI':
        attributes = computeTotalVI(configuration)
    else:
        raise ValueError(model + ' is not one of the supported models')
    configuration["attributes"] = attributes
    if get_from_config(configuration, parameters.RUN_ASYNCHRONOUSLY):
        utils.notify_backend(get_from_config(configuration, parameters.WEBHOOK), configuration)
    return configuration
# query('data/ref/source_data.h5ad', 'data/query/target_data.h5ad', 'data/model', 'data/surgery/', 'scVI')
# query(None)
