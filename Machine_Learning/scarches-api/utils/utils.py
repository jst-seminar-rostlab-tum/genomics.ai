import csv
import os
import tempfile
import scanpy
import requests


def write_latent_csv(latent, filename=tempfile.mktemp(), drop_columns=None):
    if drop_columns is None:
        drop_columns = []
    final = latent.obs.drop(columns=drop_columns)
    final["x"] = list(map(lambda p: p[0], latent.obsm["X_umap"]))
    final["y"] = list(map(lambda p: p[1], latent.obsm["X_umap"]))
    final.to_csv(filename)
    return filename


def print_csv(filename):
    with open(filename, mode='r') as file:
        for row in file.readlines():
            print(row, end='')


def save_umap_as_pdf(latent, filepath, color=None, wspace=0.6):
    if color is None:
        color = []
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
