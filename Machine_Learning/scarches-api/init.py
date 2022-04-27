from scANVI.scANVI import query
from scVI.scVI import compute_scVI
from scANVI.scANVI import compute_scANVI
import sys
import utils

def default_config():
    """
    returns the default config combined for all the models
    :return: dict containing all the default values
    """
    return {
        'model': 'scANVI',

        'reference_dataset': 'data/ref/source_data.h5ad',
        'query_dataset': 'data/query/target_data.h5ad',
        'model_path': 'data/model',
        'surgery_path': 'data/surgery',

        'pre_trained_scVI': False,
        'pre_trained_totalVI': False,
        'pre_trained_scANVI': False,

        # scANVI stuff
        'both': False,
        'compare': False,
        'predict': True,

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
        'scanvi_max_epochs': 20,
        'scvi_max_epochs': 20,
        'n_neighbors': 8,
        'max_epochs': 100,
        'unwanted_labels': ['leiden'],
        'debug': True, 
        'attributes' : None
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
    model = get_from_config(configuration, 'model')
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
    utils.notify_backend(get_from_config(configuration, 'webhook'), configuration)
# query('data/ref/source_data.h5ad', 'data/query/target_data.h5ad', 'data/model', 'data/surgery/', 'scVI')
# query(None)
