import os
import tempfile
from utils import parameters
import requests
import boto3
from aiohttp import ClientError
import scanpy
from pathlib import Path
from scarches.dataset.trvae.data_handling import remove_sparsity

UNWANTED_LABELS = ['leiden', '', '_scvi_labels', '_scvi_batch']


def get_from_config(configuration, key):
    if key in configuration:
        return configuration[key]
    return None


def to_drop(adata_obs):
    drop_list = []
    print(adata_obs)
    for attr in UNWANTED_LABELS:
        if attr in adata_obs:
            drop_list.append(attr)
    print(drop_list)
    return drop_list


def write_latent_csv(latent, key=None, filename=tempfile.mktemp(), drop_columns=None):
    """
    stores a given latent in a file, and if a key is given also in an s3 bucket
    :param latent: data to be saved
    :param key: s3 key
    :param filename: local filename
    :param drop_columns: not needed columns
    :return:
    """

    drop_columns = to_drop(latent.obs_keys())
    # if drop_columns is None:
    #     drop_columns = []
    final = latent.obs.drop(columns=drop_columns)
    final["x"] = list(map(lambda p: p[0], latent.obsm["X_umap"]))
    final["y"] = list(map(lambda p: p[1], latent.obsm["X_umap"]))
    final.to_csv(filename)
    if key is not None:
        store_file_in_s3(filename, key)
    return filename


def write_full_adata_to_csv(model, source_adata, target_adata, key=None, filename=tempfile.mktemp(), drop_columns=None,
                            cell_type_key='', condition_key='', neighbors=8):
    adata_full = source_adata.concatenate(target_adata)
    return write_adata_to_csv(model, adata_full, key, filename, drop_columns, cell_type_key, condition_key, neighbors)


def write_adata_to_csv(model, adata=None, key=None, filename=tempfile.mktemp(), drop_columns=None,
                       cell_type_key='cell_type',
                       condition_key='study', neighbors=8):
    anndata = None
    if adata is None:
        anndata = scanpy.AnnData(model.get_latent_representation())
    else:
        anndata = scanpy.AnnData(model.get_latent_representation(adata=adata))
    latent = anndata
    latent.obs['cell_type'] = adata.obs[cell_type_key].tolist()
    latent.obs['batch'] = adata.obs[condition_key].tolist()
    latent.obs['type'] = adata.obs['type'].tolist()
    scanpy.pp.neighbors(latent, n_neighbors=neighbors)
    scanpy.tl.leiden(latent)
    scanpy.tl.umap(latent)
    return write_latent_csv(latent, key, filename)


def write_combined_csv(latent_ref, latent_query, key=None, filename=tempfile.mktemp(), drop_columns=None):
    """
    stores a given latent in a file, and if a key is given also in an s3 bucket
    :param latent_ref: reference_data to be saved
    :param latent_query: query data to be saved
    :param key: s3 key
    :param filename: local filename
    :param drop_columns: not needed columns
    :return:
    """
    if drop_columns is None:
        drop_columns = []
    query = latent_query.obs.drop(columns=drop_columns)
    query["x"] = list(map(lambda p: float(p[0]), latent_query.obsm["X_umap"]))
    query["y"] = list(map(lambda p: float(p[1]), latent_query.obsm["X_umap"]))
    query["is_reference"] = ['No'] * len(latent_query.obsm["X_umap"])
    query.to_csv(filename)
    reference = latent_ref.obs.drop(columns=drop_columns)
    reference["x"] = list(map(lambda p: float(p[0]), latent_ref.obsm["X_umap"]))
    reference["y"] = list(map(lambda p: float(p[1]), latent_ref.obsm["X_umap"]))
    reference["is_reference"] = ['Yes'] * len(latent_ref.obsm["X_umap"])
    reference.to_csv(filename, header=False, mode='a')

    if key is not None:
        store_file_in_s3(filename, key)


def print_csv(filename):
    with open(filename, mode='r') as file:
        for row in file.readlines():
            print(row, end='')


def save_umap_as_pdf(latent, filepath, color=None, wspace=0.6):
    if color is None:
        color = []
    Path(os.path.dirname(filepath)).mkdir(parents=True, exist_ok=True)
    scanpy.pl.umap(latent,
                   color=color,
                   frameon=False,
                   wspace=wspace,
                   show=False,
                   save=True
                   )
    os.rename('figures/umap.pdf', filepath)


def notify_backend(endpoint, payload):
    """
    makes a post request to an endpoint specified by backend to notify them about the computed results
    :param endpoint: url
    :param payload: configuration initially specified from backend, allows them to identify which result is ready
    :return:
    """
    requests.post(endpoint, data=payload)


def fetch_file_from_s3(key, path):
    """
    downloads a file identified by a given key to a given path
    :param key: key in s3
    :param path: desired path
    :return:
    """
    client = boto3.client('s3', endpoint_url=os.getenv('AWS_ENDPOINT'),
                          aws_access_key_id=os.getenv('AWS_ACCESS_KEY'),
                          aws_secret_access_key=os.getenv('AWS_SECRET_KEY'))
    client.download_file(os.getenv('AWS_BUCKET'), key, path)


def store_file_in_s3(path, key):
    """
    stores a file in the given path in an s3 bucket
    :param path: path in our filesystem
    :param key: key to where to store the file in s3
    :return: returns ContentLength if successfully uploaded, 0 otherwise
    """
    try:
        bucket = os.getenv('AWS_BUCKET')
        client = boto3.client('s3', endpoint_url=os.getenv('AWS_ENDPOINT'),
                              aws_access_key_id=os.getenv('AWS_ACCESS_KEY'),
                              aws_secret_access_key=os.getenv('AWS_SECRET_KEY'))
        client.upload_file(path, bucket, key)
        response = client.head_object(Bucket=bucket, Key=key)
        return response['ContentLength']
    except ClientError as e:
        print(e)
    return 0


def delete_file(file):
    """
    deletes a file is found, otherwise does nothing
    :param file: file to delete
    :return:
    """
    if os.path.isfile(file):
        os.remove(file)


def read_h5ad_file_from_s3(key):
    """
    downloads an .h5ad file from s3, reads the data and deletes the file
    :param key:
    :return:
    """
    if key is None or len(key) == 0:
        return None
    filename = tempfile.mktemp(suffix=".h5ad")
    fetch_file_from_s3(key, filename)
    data = scanpy.read(filename)
    delete_file(filename)
    return data


def check_model_atlas_compatibility(model, atlas):
    compatible_atlases = []
    if model == 'scVI':
        compatible_atlases = ['pancreas', 'heart', 'human-lung', 'retina', 'fetal-immune']
    if model == 'scANVI':
        compatible_atlases = ['pancreas', 'heart', 'human-lung', 'retina', 'fetal-immune']
    if model == 'totalVI':
        compatible_atlases = ['pmbc', 'bone-marrow']
    return atlas in compatible_atlases


def pre_process_data(configuration):
    source_adata = read_h5ad_file_from_s3(get_from_config(configuration, parameters.REFERENCE_DATA_PATH))
    target_adata = read_h5ad_file_from_s3(get_from_config(configuration, parameters.QUERY_DATA_PATH))
    source_adata.obs["type"] = "reference"
    target_adata.obs["type"] = "query"
    source_adata.raw = source_adata
    try:
        source_adata = remove_sparsity(source_adata)
    except Exception as e:
        pass
    try:
        target_adata = remove_sparsity(target_adata)
    except Exception as e:
        pass
    return source_adata, target_adata


def translate_atlas_to_directory(configuration):
    atlas = get_from_config(configuration, 'atlas')
    if atlas == 'Pancreas':
        return 'pancreas'
    elif atlas == 'PBMC':
        return 'pbmc'
    elif atlas == 'Heart cell atlas':
        return 'heart'
    elif atlas == 'Human lung cell atlas':
        return 'human-lung'
    elif atlas == 'Bone marrow':
        return 'bone-marrow'
    elif atlas == 'Retina atlas':
        return 'retina'
    elif atlas == 'Fetal immune atlas':
        return 'fetal-immune'
