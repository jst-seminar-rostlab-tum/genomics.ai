# scarches API

This directory contains the code for our scarches REST-API that allows us to
provide a unified endpoint for the different scarches models to backend.

## Usage
The Dockerfile creates an image that runs a webserver which (currenty, still under
development and prone to change) *two* endpoints. `liveness` to check if the webserver is 
responsive and `query` which is the actual endpoint used to compute the `.tsv` files.

### `/liveness`
Allows `GET` requests and simply returns `up` if the webserver is working.

### `/query`
The query endpoint allows `POST` requests and takes a **valid** json which contains 
all the information for us to compute the model and store its results. We, as the ML and 
Visualisations team don't want to worry about the storage infrastructure of the backend, thus we 
expect paths that we can access to store the output, an example request body looks like this
```
{
    "surgery_path": "/dev/null",
    "model_path": "/dev/urandom",
    "generated_output_base_path": "/job/results/",
    "reference_dataset_path": "/dev/*.h5ad",
    "query_dataset_path": "/dev/*.h5ad",
    "model": "scANVI|scVI|totalVI",
    "webhook": "https://finished.computing"
    "debug": true
}
```
The given configuration is then merged with our default configuration and given to the 
models. The models then store the computed results under `generated_output_path` and 
our REST API returns the used configuration (which is previously merged with the default configuration) and computes the results asynchronously. After completing the 
query function and generating the results in the specified directories, the API will make a `POST` call 
to an endpoint specified in the original configuration as `webhook`. This request will contain 
the configuration that was used to calculate the results and allows the backend team to identify the results.

### Deployment
The Dockerfile in this directory can be used to build a deployable API image containing
the REST-API and the scarches code. It is configurable with several environment variables if needed.
- WORKERS sets the number of gunicorn workers in the container, default is 3.
- THREADS sets the number of threads, default is 8.
- PORT sets the port the container should listen on, defaults to 8080.
- AWS_BUCKET storage bucket in s3
- AWS_ENDPOINT aws endpoint
- AWS_ACCESS_KEY access key
- AWS_SECRET_KEY secret key