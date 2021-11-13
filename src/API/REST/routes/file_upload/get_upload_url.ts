import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";
import {projectModel} from "../../../../database/models/project";
import {UploadPartRequest} from "aws-sdk/clients/s3";
import s3 from "./s3";
import express from "express";
import bodyParser from "body-parser";

export default function upload_get_upload_url_route() {
    let router = express.Router();
    router.use(bodyParser.json());
    router.get('/file_upload/get_upload_url', check_auth(), async (req: ExtRequest, res) => {
        try {
            if (process.env.S3_BUCKET_NAME && req.user_id && req.query.uploadId) {
                let project = await projectModel.findOne({
                    owner: req.user_id,
                    uploadId: String(req.query.uploadId)
                }).exec();
                if (project) {
                    let params: UploadPartRequest = {
                        Bucket: process.env.S3_BUCKET_NAME,
                        Key: String(req.query.fileName),
                        PartNumber: Number(req.query.partNumber),
                        UploadId: String(req.query.uploadId)
                    }
                    let presignedUrl = await s3.getSignedUrlPromise('uploadPart', params);
                    res.status(200).send({presignedUrl});
                } else res.status(400).send("Upload was not started");
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    })
    return router;
}