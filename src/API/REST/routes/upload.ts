import express from "express";
import BluebirdPromise from "bluebird";
import AWS from "aws-sdk";
import bodyParser from "body-parser";
import check_auth from "../middleware/check_auth";

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

app.get('/start-upload', check_auth(), async (req, res, next) => {
    try {
        let params = {
            Bucket: BUCKET_NAME,
            Key: req.query.fileName,
            ContentType: req.query.fileType
        }
        let createUploadPromised = BluebirdPromise.promisify(s3.createMultipartUpload.bind(s3));
        let uploadData = await createUploadPromised(params);
        res.send({uploadId: uploadData.UploadId});
    } catch (err) {
        console.log(err);
    }
})

app.get('/get-upload-url', check_auth(), async (req, res, next) => {
    try {
        let params = {
            Bucket: BUCKET_NAME,
            Key: req.query.fileName,
            PartNumber: req.query.partNumber,
            UploadId: req.query.uploadId
        }
        console.log(params)
        let uploadPartPromised = BluebirdPromise.promisify(s3.getSignedUrl.bind(s3));
        let presignedUrl = await uploadPartPromised('uploadPart', params);
        res.send({presignedUrl});
    } catch (err) {
        console.log(err);
    }
})

app.post('/complete-upload', check_auth(), async (req, res, next) => {
    try {
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
    } catch (err) {
        console.log(err);
    }
})

app.listen(port, () => console.log(`Example app listening on port ${port}!`))