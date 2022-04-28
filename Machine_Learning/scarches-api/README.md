# scarches API

This directory contains the code for our scarches REST-API that allows us to
provide a unified endpoint for the different scarches models to backend.

## Usage

The Dockerfile creates an image that runs a webserver which (currenty, still under
development and prone to change) *two* endpoints. `liveness` to check if the webserver is
responsive and `query` which is the actual endpoint used to compute the `.tsv` files.

### `/query`

The query endpoint allows `POST` requests and takes a **valid** json which contains
all the information for us to compute the model and store its results. We, as the ML and
Visualisations team don't want to worry about the storage infrastructure of the backend, thus we
expect s3 keys that we can access to store the output. The used parameter keys that can be set in the
configuration are specified in [parameter.py](./utils/parameters.py) with a short explanation what
the an example request body looks like this

```
{
    "model": "scANVI|scVI|totalVI",
    "output_path": "[s3 key where we should store the generated .tsv]",
    "model_path": "[s3 key where we should store the resulting model.pt]",
    "ref_path": "[s3 key where we can retrieve a pretrained model]",
    "reference_data": "[s3 key to the .h5ad file used as a reference]",
    "query_data": "[s3 key to the .h5ad file used for the query]",
    "webhook": "[url where we should make a request when we complete the query]"
    "debug": true,
    "async" false
}
```

The given configuration is then merged with our default configuration and given to the
models, specified [here](./init.py). The models then store the computed results under the s3 key given by `output_path`
and our REST API returns the used configuration (which is previously merged with the default configuration) and computes
the results asynchronously. If the config contains the key `async` the computation is done
asynchronously. After completing the query function, generating the results, and storing them in the s3 buckets, the API
will make a `POST` call to an endpoint specified in the original configuration as `webhook`. This request will contain
the configuration that was used to calculate the results and allows the backend team to identify the results. If the
configuration
is done synchronously, the used configuration will be returned after the computation as a respsonse
to the original request.

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