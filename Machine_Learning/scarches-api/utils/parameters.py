"""
This file contains all the keys in the config that are changeable through the API and used throughout the computation
Each constant should be described to know what effects it has on the model.
This file does not contain the default values of the constants, this is just used so that we
can rapidly change the naming of the keys in the configuration.
"""
# sets the model
MODEL = 'model'
# sets the condition key
CONDITION_KEY = 'condition_key'
# set the cell type key
CELL_TYPE_KEY = 'cell_type_key'
# sets the target conditions
TARGET_CONDITIONS = 'target_conditions'
# sets the number of layers used
NUMBER_OF_LAYERS = 'n_layers'
# encode_covariates
ENCODE_COVARIATES = 'encode_covariates'
# deeply_inject_covariates
DEEPLY_INJECT_COVARIATES = 'deeply_inject_covariates'
# use_layer_norm
USE_LAYER_NORM = 'use_layer_norm'
# use_batch_norm
USE_BATCH_NORM = 'use_batch_norm'
# sets the identifier for the unlabeled keys
UNLABELED_KEY = 'unlabeled_key'
# sets the number of max epochs for scanvi
SCANVI_MAX_EPOCHS = 'scanvi_max_epochs'
# sets the number of max epochs for scvi
SCVI_MAX_EPOCHS = 'scvi_max_epochs'
# sets the number of neighbors
NUMBER_OF_NEIGHBORS = 'n_neighbors'
# sets the maximum number of epochs
MAX_EPOCHS = 'n_neighbors'
# sets the labeles that are not of interest and should be removed
UNWANTED_LABELS = 'unwanted_labels'
# sets if debug messages should be printed
DEBUG = 'debug'
# sets the attributes
ATTRIBUTES = 'attributes'
# compare reference and query
SCANVI_COMPARE_REFERENCE_AND_QUERY = 'scanvi_compare_reference_and_query'
# compare observed and predicted celltypes
SCANVI_COMPARE_OBSERVED_AND_PREDICTED_CELLTYPES = 'scanvi_compare_observed_and_predicted_celltypes'
# predict celltypes
SCANVI_PREDICT_CELLTYPES = 'predict'
# sets the key for the webhook to call after the computation
WEBHOOK = 'webhook'
# sets the path/s3 key of the pretrained model
PRETRAINED_MODEL_PATH = 'ref_path'
# sets the path/s3 key of the model we can store after training it
RESULTING_MODEL_PATH = 'model_path'
# sets if we can use a pretrained scVI model
USE_PRETRAINED_SCVI_MODEL = 'pre_trained_scVI'
# sets if we can use a pretrained scANVI model
USE_PRETRAINED_SCANVI_MODEL = 'pre_trained_scANVI'
# sets if we can use a pretrained totalVI model
USE_PRETRAINED_TOTALVI_MODEL = 'pre_trained_totalVI'
# sets the path/s3 key of the reference_data
REFERENCE_DATA_PATH = 'reference_data'
# sets the path/s3 key of the query data
QUERY_DATA_PATH = 'query_data'
# sets the path/s3 key of the output tsv file
OUTPUT_PATH = 'output_path'
# sets if the computation should run asynchronously
RUN_ASYNCHRONOUSLY = 'async'
