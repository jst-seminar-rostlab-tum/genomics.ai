import express from "express";
import BluebirdPromise from "bluebird";
import AWS, {S3} from "aws-sdk";
import bodyParser from "body-parser";
import check_auth from "../middleware/check_auth";
import {IProject, projectModel} from "../../../database/models/project";
import {ExtRequest} from "../../../definitions/ext_request";
import {UploadPartRequest} from "aws-sdk/clients/s3";


const app = express()
app.use(bodyParser.json())

const port = 4000
const BUCKET_NAME = process.env.S3_BUCKET_NAME;
const S3_OPTIONS = {
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

app.get('/', check_auth(), (req, res, next) => {
    res.send('Hello World!');
})

app.get('/start_upload', check_auth(), async (req: ExtRequest, res, next) => {
    try {
        if (BUCKET_NAME !== undefined) {
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

app.get('/get_upload_url', check_auth(), async (req: ExtRequest, res, next) => {
    try {
        if (BUCKET_NAME !== undefined && req.query.uploadId !== undefined) {
                let project = await projectModel.findOne({
                    owner: req.user_id,
                    uploadId: req.query.uploadId
                }).exec();
            let params: UploadPartRequest = {
                Bucket: BUCKET_NAME,
                Key: req.query.fileName,
                PartNumber: Number(req.query.partNumber),
                UploadId: req.query.uploadId
            }
            console.log(params);
            let presignedUrl = await s3.getSignedUrlPromise('uploadPart', params);
            res.send({presignedUrl});

        } else res.status(500).send("S3-BucketName is not set");
    } catch (err) {
        console.log(err);
    }
})

app.post('/complete_upload', check_auth(), async (req, res, next) => {
    try {
        if (BUCKET_NAME !== undefined) {
            console.log(req.body, ': body')
            let params = {
                Bucket: BUCKET_NAME,
                Key: req.body.params.fileName,
                MultipartUpload: {
                    Parts: req.body.params.parts
                },
                UploadId: req.body.params.uploadId
            }
            console.log(params)
            let completeUploadPromised = BluebirdPromise.promisify(s3.completeMultipartUpload.bind(s3));
            let data = await completeUploadPromised(params);
            res.send({data});
        } else res.status(500).send("S3-BucketName is not set");
    } catch (err) {
        console.log(err);
    }
})

app.listen(port, () => console.log(`Example app listening on port ${port}!`))