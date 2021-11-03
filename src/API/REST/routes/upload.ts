import express from "express";
import AWS, {AWSError, S3} from "aws-sdk";
import bodyParser from "body-parser";
import check_auth from "../middleware/check_auth";
import {IProject, projectModel} from "../../../database/models/project";
import {ExtRequest} from "../../../definitions/ext_request";
import {CompleteMultipartUploadOutput, CompleteMultipartUploadRequest, UploadPartRequest} from "aws-sdk/clients/s3";


const app = express()
app.use(bodyParser.json())

const port = 4000
const BUCKET_NAME = process.env.S3_BUCKET_NAME;
const S3_OPTIONS: S3.Types.ClientConfiguration = {
    accessKeyId: process.env.S3_ACCESS_KEY_ID,
    secretAccessKey: process.env.S3_SECRET_ACCESS_KEY,
    endpoint: process.env.S3_ENDPOINT,
    s3ForcePathStyle: true,
    signatureVersion: 'v4'
};
const s3 = new AWS.S3(S3_OPTIONS);

app.use(function (req, res, next) {
    res.header("Access-Control-Allow-Origin", "*");
    res.header("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept");
    next();
});


app.get('/start_upload', check_auth(), async (req: ExtRequest, res) => {
    try {
        if (BUCKET_NAME  && req.query.fileName) {
            let project: IProject = await projectModel.create({
                owner: req.user_id,
                fileName: req.query.fileName,
                uploadDate: new Date(),
                status: "UPLOAD_PENDING"
            });
            let params: S3.CreateMultipartUploadRequest = {
                Bucket: BUCKET_NAME,
                Key: project._id
            }
            s3.createMultipartUpload(params, (err, uploadData) => {
                if (err)
                    console.error(err, err.stack || "Error when requesting uploadId");
                else
                    res.send({uploadId: uploadData.UploadId});
            });
        } else res.status(500).send("S3-BucketName is not set");
    } catch (err) {
        console.log(err);
    }
})

app.get('/get_upload_url', check_auth(), async (req: ExtRequest, res) => {
    try {
        if (BUCKET_NAME  && req.query.uploadId ) {
            let project = await projectModel.findOne({
                owner: req.user_id,
                uploadId: String(req.query.uploadId)
            }).exec();
            if (project) {
                let params: UploadPartRequest = {
                    Bucket: BUCKET_NAME,
                    Key: project._id,
                    PartNumber: Number(req.query.partNumber),
                    UploadId: String(req.query.uploadId)
                }
                console.log(params);
                let presignedUrl = await s3.getSignedUrlPromise('uploadPart', params);
                res.send({presignedUrl});
            } else res.status(400).send("Upload was not started");
        } else res.status(500).send("S3-BucketName is not set");
    } catch (err) {
        console.log(err);
    }
})

app.post('/complete_upload', check_auth(), async (req: ExtRequest, res) => {
    try {
        if (BUCKET_NAME ) {
            console.log(req.body, ': body')
            let params: CompleteMultipartUploadRequest = {
                Bucket: BUCKET_NAME,
                Key: req.body.params.fileName,
                MultipartUpload: {
                    Parts: req.body.params.parts
                },
                UploadId: req.body.params.uploadId
            }
            console.log(params)
            s3.completeMultipartUpload(params, (err: AWSError, data: CompleteMultipartUploadOutput) => {
                if (err) console.error(err, err.stack || "Error when completing multipart upload");
                if (data.Key  && data.Bucket ) {
                    let request: S3.HeadObjectRequest = {Key: data.Key, Bucket: data.Bucket};
                    s3.headObject(request)
                        .promise()
                        .then(result => projectModel.updateOne({uploadId: params.UploadId},
                            {fileSize: result.ContentLength}));
                }
                if (data.Location)
                    projectModel.updateOne({uploadId: params.UploadId},
                        {location: data.Location});
                res.send({data});
            });
        } else res.status(500).send("S3-BucketName is not set");
    } catch (err) {
        console.log(err);
    }
})

app.listen(port, () => console.log(`Example app listening on port ${port}!`))