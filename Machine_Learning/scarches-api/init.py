from scANVI.scANVI import query
from scVI.scVI import compute_scVI
from scANVI.scANVI import compute_scANVI
import sys
from utils import utils, parameters


def default_config():
    """
    returns the default config combined for all the models
    :return: dict containing all the default values
    """
    return {
        parameters.MODEL: 'scANVI',

        'reference_dataset': 'data/ref/source_data.h5ad',
        'query_dataset': 'data/query/target_data.h5ad',
        'model_path': 'data/model',
        'surgery_path': 'data/surgery',

        'pre_trained_scVI': False,
        'pre_trained_totalVI': False,
        'pre_trained_scANVI': False,

        # scANVI stuff
        parameters.SCANVI_COMPARE_REFERENCE_AND_QUERY: False,
        parameters.SCANVI_COMPARE_OBSERVED_AND_PREDICTED_CELLTYPES: False,
        'predict': True,

        parameters.CONDITION_KEY: 'study',
        parameters.CELL_TYPE_KEY: 'cell_type',
        parameters.TARGET_CONDITIONS: ['Pancreas CelSeq2', 'Pancreas SS2'],
        'ref_path': './ref_model/',
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
        parameters.DEBUG: True,
        parameters.ATTRIBUTES: None
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
        # compute_totalVI()
        print("not supported yet")
    else:
        raise ValueError(model + ' is not one of the supported models')
    configuration["attributes"] = attributes
    utils.notify_backend(get_from_config(configuration, parameters.WEBHOOK), configuration)
# query('data/ref/source_data.h5ad', 'data/query/target_data.h5ad', 'data/model', 'data/surgery/', 'scVI')
# query(None)
