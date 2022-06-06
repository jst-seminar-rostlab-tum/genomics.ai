# scarches API

ScArches is a novel deep learning model that enables mapping query to reference datasets. The model allows the user to construct single or multi-modal (CITE-seq) references as well as classifying unlabelled query cells.
Currently, mapping of query to reference datasets is offered on GeneCruncher.

This directory [scarches-api](../scarches-api) contains the code for our scarches REST-API that allows us to
provide a unified endpoint for the different scarches models to backend
and compute a `query.csv` file that can be downloaded from the backend s3 bucket and parse for the visualization.

***
## Models and workflow
Our code support the following models: scVI, scANVI and totalVI.

### scVI
[scVI](./scVI/scVI.py) is an unsupervised model that maps query to reference atlases. It does not require cell type labels for the mapping.

The workflow is: 
* set up modules and figure parameters
* [download](./utils/utils.py) reference and query dataset from s3
* preprocess reference and query dataset
* create scVI model (or use pretrained model) and train it on reference dataset
* create anndata file of latent representation and compute UMAP
* train the model on query dataset and save the result
* [write](./utils/utils.py) the result into `query.csv` file and [upload](./utils/utils.py) to s3

### scANVI
[scANVI](./scANVI/scANVI.py) is a semi-supervised variant of scVI designed to leverage any available cell state annotations. Compared to unsupervised models, this model will perform better integration if cell type labels are partially available in the query.

The workflow is:
* set up modules and figure parameters
* [download](./utils/utils.py) reference and query dataset from s3
* preprocess reference and query dataset
* create scANVI model (or use pretrained model) and train it on reference dataset
* create anndata file of latent representation and compute UMAP 
* predict unlabeled cell type of latent and compute the accuracy of the learned classifier
* train the model on query dataset and save the result
* [write](./utils/utils.py) the result into `query.csv` file and [upload](./utils/utils.py) to s3

### totalVI
[totalVI](./totalVI/totalVI.py) is a multi-modal model that can be used to map to multi-modal reference atlases. It can be used for the imputation of proteins in the query dataset.

The workflow is:
* set up modules and figure parameters
* [download](./utils/utils.py) reference and query dataset from s3
* preprocess reference and query dataset
* set up anndata and create totalVI model or use pretrained model 
* train model on reference dataset
* save latent representation and visualize RNA data
* perform surgery on reference model and train on query dataset without protein data
* impute protein data for the query dataset and visualize result
* get latent representation of reference + query dataset and compute UMAP
* [write](./utils/utils.py) the result into `query.csv` file and [upload](./utils/utils.py) to s3
***

## Atlas
Our code supports the following Atlases. Their compatible models are: 
- Pancreas (scVI, scANVI)
- PBMC (totalVI)
- Heart cell atlas (scVI, scANVI)
- Human lung cell atlas (scANVI)
- Retina atlas (scVI, scANVI)
- Fetal immune atlas (scVI, scANVI)

***
## Usage

The Dockerfile creates an image that runs a webserver which provide *two* endpoints: 
- `liveness` to check if the webserver is
responsive and 
- `query` which is the actual endpoint used to compute the `.csv` files, which is used for the visualization.

### `/query`

The query endpoint allows `POST` requests and takes a **valid** json which contains
all the information for us to compute the model and store its results. On the ML and
Visualisation side we expect only s3 keys (from the backend) that we can access to store the output.
The used parameter keys that can be set in the
configuration are specified in [parameter.py](./utils/parameters.py) with a short explanation what
the parameter does and what the key should look like. 

The given configuration is then merged with our default configuration and given to the
models, specified in [init.py](./init.py). The models then store the computed results under the s3 key given by `output_path`
and our REST-API returns the used configuration (which is previously merged with the default configuration). 
If the configuration contains the key `async` the computation is done
_asynchronously_. After completing the query function, generating the results, and storing them in the s3 buckets, the API
will make a `POST` request to an endpoint specified in the original configuration as `webhook`. This request contains all information about
the configuration that was used to calculate the results and allows the backend team to identify the results. 
If the configuration is done _synchronously_, the given configuration will be returned after the computation as a response
to the original request.

#### Example request body

```
{
    "model": "scANVI|scVI|totalVI",
    "atlas": "Pancreas|PBMC|Heart cell atlas|Human lung cell atlas|Bone marrow|Retina atlas|Fetal immune atlas",
    "output_path": "[s3 key where we should store the generated .csv]",
    "reference_data": "[s3 key to the .h5ad file used as a reference]",
    "query_data": "[s3 key to the .h5ad file used for the query]",
    "webhook": "[url where we should make a request when we complete the query]",
    "async": false
}
```

#### Example with scVI
##### Atlas: Pancreas
request
```
{
    "model": "scVI",
    "atlas": "Pancreas"
    "output_path": "query.csv",
    "reference_data": "pancreas_source.h5ad",
    "query_data": "pancreas_query.h5ad"
}
```
response
```
{
    "async": false,
    "attributes": null,
    "cell_type_key": "cell_type",
    "condition_key": "study",
    "debug": false,
    "deeply_inject_covariates": false,
    "development-mode": false,
    "encode_covariates": true,
    "max_epochs": 100, 
    "model": "scVI",
    "model_path": "model.pt",
    "n_layers": 2,
    "n_neighbors": 8,
    "output_path": "query.csv",
    "pre_trained_scANVI": true,
    "pre_trained_scVI": true,
    "pre_trained_totalVI": true,
    "predict": false,
    "query_data": "pancreas_query.h5ad",
    "ref_path": "model.pt",
    "reference_data": "pancreas_source.h5ad",
    "scanvi_do_surgery": false,
    "scanvi_max_epochs": 20,
    "scanvi_query_max_epochs": 100,
    "scvi_max_epochs": 400,
    "scvi_query_max_epochs": 200,
    "unlabeled_key": "Unknown",
    "unwanted_labels": [
        "leiden"
    ],
    "use_gpu": false,
    "use_batch_norm": "none",
    "use_layer_norm": "both",
    "webhook": ""
}
```

#### Example with scANVI
##### Atlas: Human lung cell
request
```
{
    "model": "scANVI",
    "atlas": "Human lung cell atlas"
    "output_path": "query.csv",
    "model_path": "model.pt",
    "pre_trained_scVI": false,
    "reference_data": "hcla_working.h5ad",
    "query_data": "human_lung_query_working.h5ad",
    "webhook": "",
    "ref_path": "model.pt",
    "debug": false,
    "async": false
}

```
response
```
{
    "async": false,
    "attributes": null,
    "cell_type_key": "scanvi_label",
    "condition_key": "dataset",
    "debug": false,
    "deeply_inject_covariates": false,
    "encode_covariates": true,
    "max_epochs": 100, 
    "model": "scANVI",
    "model_path": "model.pt",
    "n_layers": 2,
    "n_neighbors": 8,
    "output_path": "query.csv",
    "pre_trained_scANVI": false,
    "pre_trained_scVI": true,
    "pre_trained_totalVI": false,
    "predict": true,
    "query_data": "human_lung_query_working.h5ad",
    "ref_path": "model.pt",
    "reference_data": "hcla_working.h5ad",
    "scanvi_max_epochs": 20,
    "scanvi_query_max_epochs": 100,
    "scvi_max_epochs": 400,
    "scvi_query_max_epochs": 200,    
    "unlabeled_key": "unlabled",
    "unwanted_labels": [
        "leiden"
    ],
    "use_batch_norm": "none",
    "use_layer_norm": "both",
    "webhook": ""
}
```

#### Example with totalVI
##### Atlas: PBMC

request
```
{
    "model": "totalVI",
    "atlas": "PBMC",
    "output_path": "query.csv",
    "reference_data": "pbmc_reference.h5ad",
    "query_data": "pbmc_query.h5ad",
}
```
response
```
{
    "async": false,
    "attributes": null,
    "cell_type_key": "cell_type",
    "condition_key": "study",
    "debug": false,
    "deeply_inject_covariates": false,
    "encode_covariates": true,
    "max_epochs": 100, 
    "model": "totalVI",
    "model_path": "model.pt",
    "n_layers": 2,
    "n_neighbors": 8,
    "output_path": "query.csv",
    "pre_trained_scANVI": false,
    "pre_trained_scVI": false,
    "pre_trained_totalVI": false,
    "predict": false,
    "query_data": "pbmc_query.h5ad",
    "ref_path": "model.pt",
    "reference_data": "pbmc_reference.h5ad",
    "scanvi_max_epochs": 20,
    "scanvi_query_max_epochs": 100,
    "scvi_max_epochs": 400,
    "scvi_query_max_epochs": 200,
    "totalvi_max_epochs_1": 400,
    "totalvi_max_epochs_2": 200,
    "unlabeled_key": "Unknown",
    "unwanted_labels": [
        "leiden"
    ],
    "use_batch_norm": "none",
    "use_layer_norm": "both",
    "webhook": ""
}
```



### Deployment

The Dockerfile in this directory can be used to build a deployable API image containing
the REST-API and the scarches code with the above-mentioned models. It is configurable with several environment variables if needed.

- WORKERS sets the number of gunicorn workers in the container, default is 3.
- THREADS sets the number of threads, default is 8.
- PORT sets the port the container should listen on, defaults to 8080.
- AWS_BUCKET storage bucket in s3
- AWS_ENDPOINT aws endpoint
- AWS_ACCESS_KEY access key
- AWS_SECRET_KEY secret key