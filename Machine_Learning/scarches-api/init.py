import os
import time

startTime = time.time()

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
        parameters.ATLAS: 'Human lung cell atlas',

        parameters.REFERENCE_DATA_PATH: 'hcla.h5ad',
        parameters.QUERY_DATA_PATH: 'last_attempt.h5ad',
        parameters.RESULTING_MODEL_PATH: 'model.pt',
        parameters.OUTPUT_PATH: 'query.csv',

        parameters.USE_PRETRAINED_SCVI_MODEL: False,
        parameters.USE_PRETRAINED_TOTALVI_MODEL: False,
        parameters.USE_PRETRAINED_SCANVI_MODEL: True,
        parameters.USE_GPU: True,

        # scANVI stuff
        parameters.SCANVI_COMPARE_REFERENCE_AND_QUERY: False,
        parameters.SCANVI_COMPARE_OBSERVED_AND_PREDICTED_CELLTYPES: False,
        parameters.SCANVI_PREDICT_CELLTYPES: True,

        parameters.CONDITION_KEY: 'study',
        parameters.CELL_TYPE_KEY: 'cell_type',
        parameters.PRETRAINED_MODEL_PATH: '',
        parameters.NUMBER_OF_LAYERS: 2,
        parameters.ENCODE_COVARIATES: True,
        parameters.DEEPLY_INJECT_COVARIATES: False,
        parameters.USE_LAYER_NORM: 'both',
        parameters.USE_BATCH_NORM: 'none',
        parameters.UNLABELED_KEY: 'Unknown',
        parameters.SCANVI_MAX_EPOCHS: 20,
        parameters.SCANVI_MAX_EPOCHS_QUERY: 100,
        parameters.SCVI_MAX_EPOCHS: 400,
        parameters.SCVI_QUERY_MAX_EPOCHS: 200,
        parameters.NUMBER_OF_NEIGHBORS: 8,
        parameters.MAX_EPOCHS: 100,
        parameters.UNWANTED_LABELS: ['leiden'],
        parameters.DEBUG: False,
        parameters.RUN_ASYNCHRONOUSLY: False,
        parameters.ATTRIBUTES: None,

        # totalVI stuff
        parameters.TOTALVI_MAX_EPOCHS_1: 400,
        parameters.TOTALVI_MAX_EPOCHS_2: 200,

        parameters.SCANVI_DO_SURGERY: False,
        parameters.DEV_DEBUG: False,
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
    start_time = time.time()
    configuration = merge_configs(user_config)
    #Sets the correct condition and cell_type key
    configuration = utils.set_keys(configuration)
    
    model = utils.get_from_config(configuration, parameters.MODEL)
    configuration['atlas'] = utils.translate_atlas_to_directory(configuration)
    if model == 'scVI':
        attributes = compute_scVI(configuration)
    elif model == 'scANVI':
        attributes = compute_scANVI(configuration)
    elif model == 'totalVI':
        attributes = computeTotalVI(configuration)
    else:
        raise ValueError(model + ' is not one of the supported models')
    configuration["attributes"] = attributes
    run_time = (time.time() - start_time)
    print('completed query in ' + str(run_time) + 's and stored it in: ' + get_from_config(configuration,
                                                                                           parameters.OUTPUT_PATH))
    if get_from_config(configuration, parameters.WEBHOOK) is not None and len(
            get_from_config(configuration, parameters.WEBHOOK)) > 0:
        utils.notify_backend(get_from_config(configuration, parameters.WEBHOOK), configuration)
    return configuration


# query('data/ref/source_data.h5ad', 'data/query/target_data.h5ad', 'data/model', 'data/surgery/', 'scVI')

if __name__ == "__main__":
    os.environ["AWS_BUCKET"] = 'minio-bucket'
    os.environ['AWS_ENDPOINT'] = 'http://127.0.0.1:9000'
    os.environ['AWS_ACCESS_KEY'] = 'minioadmin'
    os.environ['AWS_SECRET_KEY'] = 'minioadmin'

    query({parameters.DEV_DEBUG: True})

