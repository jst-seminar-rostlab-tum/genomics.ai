from distutils.command.config import config
from scANVI.scANVI import query
from scVI.scVI import compute_scVI
from scANVI.scANVI import compute_scANVI

config = {
    'model' : 'scANVI',

    'reference_dataset' : 'data/ref/source_data.h5ad',
    'query_dataset' : 'data/query/target_data.h5ad',
    'model_path' : 'data/model',
    'surgery_path' : 'data/surgery',

    'pre_trained_scVI' : False,
    'pre_trained_totalVI' : False,
    'pre_trained_scANVI' : False,

    # scANVI stuff
    'both' : False,
    'compare' : False,
    'predict' : True,

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

    'debug' : True
}

def get_from_config(key):
    
    return config[key]

#def query(reference_dataset, query_dataset, model_path, surgery_path,  model_type):
def query():
    if get_from_config('model') == 'scVI' :
        compute_scVI(config)

    elif get_from_config('model') == 'scANVI' :
        compute_scANVI(config)
    elif get_from_config('model') == 'totalVI' :
        #compute_totalVI()
        print("wrong")
    else :
        raise ValueError(get_from_config('model') + ' is not one of the supported models')

#query('data/ref/source_data.h5ad', 'data/query/target_data.h5ad', 'data/model', 'data/surgery/', 'scVI')
query()
