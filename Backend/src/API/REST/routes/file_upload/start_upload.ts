import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";
import {IProjectJob, projectJobModel} from "../../../../database/models/projectJob";
import {S3} from "aws-sdk";
import s3 from "../../../../util/s3";
import express from "express";

export default function upload_start_upload_route() {
    let router = express.Router();
    router.post('/file_upload/start_upload', check_auth(), async (req: ExtRequest, res) => {
        let {fileName} = req.body;
        if(!fileName)
            return res.status(400).send("Missing fileName parameter.");

        try {
            if (process.env.S3_BUCKET_NAME && req.user_id) {
                let project: IProjectJob = await projectJobModel.create({
                    owner: req.user_id,
                    fileName: String(fileName),
                    uploadDate: new Date(),
                    status: "UPLOAD_PENDING"
                });
                let params: S3.CreateMultipartUploadRequest = {
                    Bucket: process.env.S3_BUCKET_NAME,
                    Key: String(project._id)
                }
                s3.createMultipartUpload(params, (err, uploadData) => {
                    if (err) {
                        console.error(err, err.stack || "Error when requesting uploadId");
                        res.status(500).send(err);
                    } else {
                        projectJobModel.updateOne(
                            {_id: project._id},
                            {uploadId: uploadData.UploadId}
                        ).exec();
                        res.status(200).send({uploadId: uploadData.UploadId});
                    }
                });
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    })
    return router;
}