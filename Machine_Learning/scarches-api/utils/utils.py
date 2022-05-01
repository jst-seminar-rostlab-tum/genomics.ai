import os
import tempfile
import scanpy
import requests
import boto3
from aiohttp import ClientError
import sys
from pathlib import Path
import scanpy
import fileinput
from pathlib import Path


def write_latent_csv(latent, key=None, filename=tempfile.mktemp(), drop_columns=None):
    """
    stores a given latent in a file, and if a key is given also in an s3 bucket
    :param latent: data to be saved
    :param key: s3 key
    :param filename: local filename
    :param drop_columns: not needed columns
    :return:
    """
    if drop_columns is None:
        drop_columns = []
    final = latent.obs.drop(columns=drop_columns)
    final["x"] = list(map(lambda p: p[0], latent.obsm["X_umap"]))
    final["y"] = list(map(lambda p: p[1], latent.obsm["X_umap"]))
    final.to_csv(filename)
    if key is not None:
        store_file_in_s3(filename, key)
    return filename


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
    filename = tempfile.mktemp(suffix=".h5ad")
    fetch_file_from_s3(key, filename)
    data = scanpy.read(filename)
    delete_file(filename)
    return data
