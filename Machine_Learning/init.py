from distutils.command.config import config
from scVI.scVI import compute_scVI

config = {
    'condition_key': 'study',
    'cell_type_key': 'cell_type',
    'target_conditions': ['Pancreas CelSeq2', 'Pancreas SS2'],
    'n_layers': 2,
    'encode_covariates': True,
    'deeply_inject_covariates': False,
    'use_layer_norm': 'both',
    'use_batch_norm': 'none',
    'unlabeled_key': 'Unknown',
    'vae_max_epochs': 20,
    'n_neighbors': 8,
    'max_epochs': 100,
    'unwanted_labels': ['leiden'],
    'model_name' : 'model.pt' # To see if the model already exists
}

def query(reference_dataset, query_dataset, model_path, surgery_path,  model_type):
    if model_type == 'scVI' :
        compute_scVI(config, reference_dataset, query_dataset, surgery_path, model_path, True)
    elif model_type == 'scanVI' :
        #compute_scanVI()
        print("wrong")
    elif model_type == 'totalVI' :
        #compute_totalVI()
        print("wrong")
    else :
        raise ValueError(model_type + ' is not one of the supported models')

query('data/ref/source_data.h5ad', 'data/query/target_data.h5ad', 'data/model', 'data/surgery/', 'scVI')
