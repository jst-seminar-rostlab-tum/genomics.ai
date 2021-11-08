import express from "express";
import AWS, {AWSError, S3} from "aws-sdk";
import bodyParser from "body-parser";
import check_auth from "../middleware/check_auth";
import {IProject, projectModel} from "../../../database/models/project";
import {ExtRequest} from "../../../definitions/ext_request";
import {CompleteMultipartUploadOutput, CompleteMultipartUploadRequest, UploadPartRequest} from "aws-sdk/clients/s3";

export default function upload_route() {
    let router = express.Router();
    router.use(bodyParser.json());

    const BUCKET_NAME = process.env.S3_BUCKET_NAME;
    const S3_OPTIONS: S3.Types.ClientConfiguration = {
        accessKeyId: process.env.S3_ACCESS_KEY_ID,
        secretAccessKey: process.env.S3_SECRET_ACCESS_KEY,
        endpoint: process.env.S3_ENDPOINT,
        s3ForcePathStyle: true,
        signatureVersion: 'v4'
    };
    const s3 = new AWS.S3(S3_OPTIONS);

    router.use(function (req, res, next) {
        res.header("Access-Control-Allow-Origin", "*");
        res.setHeader('Access-Control-Allow-Methods', 'GET,POST,PUT,PATCH,DELETE');
        res.header("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept, auth");
        next();
    });


    router.get('/start_upload', check_auth(), async (req: ExtRequest, res) => {
        try {
            if (BUCKET_NAME && req.user_id && req.query.fileName) {
                let project: IProject = await projectModel.create({
                    owner: req.user_id,
                    fileName: req.query.fileName,
                    uploadDate: new Date(),
                    status: "UPLOAD_PENDING"
                });
                let params: S3.CreateMultipartUploadRequest = {
                    Bucket: BUCKET_NAME,
                    Key: String(project._id)
                }
                s3.createMultipartUpload(params, (err, uploadData) => {
                    if (err)
                        console.error(err, err.stack || "Error when requesting uploadId");
                    else {
                        projectModel.updateOne(
                            {_id: project._id},
                            {uploadId: uploadData.UploadId}
                        ).exec();
                        res.send({uploadId: uploadData.UploadId});
                    }
                });
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
        }
    })

    router.get('/get_upload_url', check_auth(), async (req: ExtRequest, res) => {
        try {
            if (BUCKET_NAME && req.user_id && req.query.uploadId) {
                let project = await projectModel.findOne({
                    owner: req.user_id,
                    uploadId: String(req.query.uploadId)
                }).exec();
                if (project) {
                    let params: UploadPartRequest = {
                        Bucket: BUCKET_NAME,
                        Key: String(project._id),
                        PartNumber: Number(req.query.partNumber),
                        UploadId: String(req.query.uploadId)
                    }
                    let presignedUrl = await s3.getSignedUrlPromise('uploadPart', params);
                    res.send({presignedUrl});
                } else res.status(400).send("Upload was not started");
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
        }
    })

    router.post('/complete_upload', check_auth(), async (req: ExtRequest, res) => {
        try {
            if (BUCKET_NAME && req.user_id) {
                let project = await
                    projectModel.findOne({
                        owner: req.user_id,
                        uploadId: String(req.body.params.uploadId)
                    }).exec();
                if (project) {
                    //console.log(req.body, ': body')
                    let params: CompleteMultipartUploadRequest = {
                        Bucket: BUCKET_NAME,
                        Key: String(project._id),
                        MultipartUpload: {
                            Parts: req.body.params.parts
                        },
                        UploadId: req.body.params.uploadId
                    }
                    //console.log(params)
                    s3.completeMultipartUpload(params, (err: AWSError, data: CompleteMultipartUploadOutput) => {
                        if (err) console.error(err, err.stack || "Error when completing multipart upload");
                        if (data.Key && data.Bucket) {
                            let request: S3.HeadObjectRequest = {Key: data.Key, Bucket: data.Bucket};
                            s3.headObject(request)
                                .promise()
                                .then(result => projectModel.updateOne({uploadId: params.UploadId},
                                    {fileSize: result.ContentLength}).exec());
                        }
                        if (data && data.Location && project) {
                            projectModel.updateOne(
                                {_id: project._id},
                                {
                                    location: data.Location
                                    , status: "UPLOAD_COMPLETE"
                                }).exec();
                        }
                        res.send({data});
                    });

                } else res.status(400).send("Project could not be found");
            } else {
                res.status(500).send("Server was not set up correctly")
            }
        } catch (err) {
            console.log(err);
        }
    })
    return router;
}