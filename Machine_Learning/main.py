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
endpoint = os.getenv('ENDPOINT')
access_key = os.getenv('ACCESS_KEY')
secret_key = os.getenv(
    'SECRET_KEY')
modelname = os.getenv('MODEL_NAME')
database_uri = os.getenv(
    'DATABASE_URI')
bucket = os.getenv('BUCKET')


def dispose(filename):
    if os.path.isfile(filename):
        print(filename + " found, will be disposed")
        os.remove(filename)


s3 = boto3.client('s3', endpoint_url=endpoint,
                  aws_access_key_id=access_key, aws_secret_access_key=secret_key)
print("starting download")
s3.download_file(bucket, modelname, modelname)
print("Download finished, loading model")
clf = joblib.load(modelname)
print("Model loaded, ready to dispose")
dispose(modelname)

db = MongoClient(database_uri).get_default_database()


@app.route("/")
@app.route("/<path:path>")
def catch_all(path):
    return 'You want path: %s, which is not yet implemented or does not exist' % path


@app.route("/run_classifier")
def classify():
    data = request.args
    for key in ["uploadId"]:
        if key not in data.keys():
            return "Key \"{}\" missing in request json data!\nPlease check again if the request is correct!".format(
                key), 400
    uploadId = data['uploadId']
    project = db.projects.find_one({"uploadId": uploadId})

    if project is None:
        message = f"There exists no project with upload_id {uploadId}"
        print(message)
        return message, 400
    if (project['status']) == "ABORTED":
        print("Project has been aborted. Terminating.")
        return
    print("Project found and not aborted")
    fileName = str(project['_id'])
    print("Starting download h5ad")
    s3.download_file(bucket, fileName, fileName + '.h5ad')
    print("Ready for prediction")
    result = predict(fileName + '.h5ad')
    uploadSize = upload(result)
    dispose(result)
    dispose(fileName + '.h5ad')

    db.projects.update_one({'uploadId': uploadId}, {
        "$set": {"status": "DONE", "resultSize": uploadSize, "resultName": result}})
    print("Classification has been computed")
    return "Classification has been computed", 200


def upload(filename):
    s3.upload_file(filename, bucket, filename)
    response = s3.head_object(Bucket=bucket, Key=filename)
    return response['ContentLength']


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
    cleanedDataset['x'] = list(
        map(lambda pair: pair[0], input.obsm['X_umap']))

    cleanedDataset['y'] = list(
        map(lambda pair: pair[1], input.obsm['X_umap']))
    resultname = 'result_' + filename.rsplit(".", 1)[0] + '.tsv'
    cleanedDataset.index.name = 'id'
    cleanedDataset.to_csv(resultname, columns=['x', 'y', 'celltype'], sep='\t')

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
