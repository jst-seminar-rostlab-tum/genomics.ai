
from __future__ import annotations
import os
from os import environ
from flask import Flask, request
import scanpy
import joblib
import pandas as pd
from pymongo import MongoClient
import boto3
requiredKeysDefault = {
    'nCount_ADT': 0, 'nFeature_ADT': 0, 'nCount_RNA': 0, 'nFeature_RNA': 0,
    'nCount_SCT': 0, 'nFeature_SCT': 0, 'Phase_G1': 0, 'Phase_G2M': 0, 'Phase_S': 1}
outputLabels = ['B', 'CD4 T', 'CD8', 'DC',
                'Mono', 'NK', 'other', 'other T']

app = Flask(__name__)


modelFile = 'model_small.joblib'
# modelFile=download_file(environ["MODEL_NAME"])
clf = joblib.load(modelFile)

db = MongoClient(environ["DATABASE_URI"]).get_default_database()


@app.route("/")
@app.route("/<path:path>")
def catch_all(path):
    return 'You want path: %s, which is not yet implemented or does not exist' % path


@app.route("/run_classifier", methods=['POST'])
def classify():
    data = request.get_json()
    for key in ["_id"]:
        if key not in data.keys():
            return "Key \"{}\" missing in request json data!\nPlease check again if the request is correct!".format(key), 400
    try:
        project_id = request.args.get('_id')
        project = db.projects.find_one({"_id": project_id})

        if project is None:
            message = f"There exists no project with _id {project_id}"
            print(message)
            return message,400
        if int(project['status']) == 3:
            print("Project has been aborted. Terminating.")
            return

        filename = download_file(project.location)
        result = predict(filename)
        upload(result)
        if os.path.isfile(result):
            print("file found")
            os.remove(result)
        return "Classification has been computed", 200
    except:
        return f"Please provide a valid bucket ID!"


def upload(filename):
    s3 = boto3.client('s3', endpoint_url=os.environ['ENDPOINT'],
                      aws_access_key_id=os.environ['ACCESS_KEY'], aws_secret_access_key=os.environ['secret_key'])
    s3.upload_file(filename, os.environ['bucket'], filename)


def predict(filename):
    input = scanpy.read_h5ad(filename, backed='r+')
    cleanedDataset = input.obs
    if not cleanedDataset.empty:
        cleanedDataset = pd.get_dummies(cleanedDataset)
    for key in cleanedDataset:
        if key not in requiredKeysDefault.keys():
            cleanedDataset.drop(key, inplace=True, axis=1)
    for key, default in requiredKeysDefault.items():
        if key not in cleanedDataset.keys():
            cleanedDataset[key] = default
    cleanedDataset = cleanedDataset[requiredKeysDefault.keys()]
    y_predict = clf.predict(cleanedDataset)
    output = pd.DataFrame(
        data=y_predict, index=cleanedDataset.index, columns=outputLabels)
    output = output.idxmax(axis=1)
    if "X_umap" not in input.obsm.keys():
        scanpy.pp.normalize_total(input)
        scanpy.pp.log1p(input)
        scanpy.pp.pca(input)
        scanpy.pp.neighbors(input)
        scanpy.tl.umap(input)

    cleanedDataset['celltype'] = output
    cleanedDataset['umap_0'] = list(
        map(lambda pair: pair[0], input.obsm['X_umap']))

    cleanedDataset['umap_1'] = list(
        map(lambda pair: pair[1], input.obsm['X_umap']))
    resultname = 'result_'+filename.rsplit(".", 1)[0]+'.tsv'
    cleanedDataset.to_csv(resultname, columns=[
                          'celltype', 'umap_0', 'umap_1'], sep='\t')

    return resultname


def download_file(url):
    print("Begin Download")
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                f.write(chunk)
    return local_filename


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
    filename = download_file(
        "https://fastdl.mongodb.org/windows/mongodb-windows-x86_64-5.0.5-signed.msi")
    print(filename)
