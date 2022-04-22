import csv
import os
import tempfile
import scanpy

def write_latent_csv(latent, filename=tempfile.mktemp(), drop_columns=[]):
    final = latent.obs.drop(columns=drop_columns)
    final["x"] = list(map(lambda p: p[0], latent.obsm["X_umap"]))
    final["y"] = list(map(lambda p: p[1], latent.obsm["X_umap"]))
    final.to_csv(filename)
    return filename


def print_csv(filename):
    with open(filename, mode='r') as file:
        for row in file.readlines():
            print(row, end='')

def save_umap_as_pdf(latent, filepath, color=[], wspace=0.6):
    scanpy.pl.umap(latent,
                   color=color,
                   frameon=False,
                   wspace=wspace,
                   show=False,
                   save=True
                   )
    os.rename('figures/umap.pdf', filepath)