import express from "express";
import {CompleteMultipartUploadOutput, CompleteMultipartUploadRequest} from "aws-sdk/clients/s3";
import {AWSError, S3} from "aws-sdk";

import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";
import {ProjectJobStatus, projectModel} from "../../../../database/models/project";
import s3 from "../../../../util/s3";
import {GoogleAuth} from "google-auth-library";

export default function upload_complete_upload_route() {
    let router = express.Router();
    router.post('/file_upload/complete_upload', check_auth(), async (req: ExtRequest, res) => {
        let {parts, uploadId} = req.body;
        if(!parts)
            return res.status(400).send("Missing parts parameter.");
        if(!uploadId)
            return res.status(400).send("Missing uploadId parameter.");

        try { // todo: all routes should be try-caught by default
            if (!process.env.S3_BUCKET_NAME || !req.user_id)
                return res.status(500).send("Server was not set up correctly");

            let project = await projectModel.findOne({
                    owner: req.user_id,
                    uploadId: String(uploadId)
                }).exec();
            if(!project)
                return res.status(400).send("Project could not be found");

            let params: CompleteMultipartUploadRequest = {
                Bucket: process.env.S3_BUCKET_NAME,
                Key: String(project._id),
                MultipartUpload: { Parts: parts },
                UploadId: String(uploadId)
            }

            s3.completeMultipartUpload(params, (err: AWSError, data: CompleteMultipartUploadOutput) => {
                if (err) {
                    console.error(err, err.stack || "Error when completing multipart upload");
                    return res.status(500).send(err);
                }
                if (!data || !data.Key || !data.Bucket || !data.Location)
                    return res.status(500).send("Error getting Multipart-Upload object data");

                let request: S3.HeadObjectRequest = {Key: data.Key, Bucket: data.Bucket};
                s3.headObject(request)
                    .promise()
                    .then(result => projectModel.updateOne(
                        {uploadId: params.UploadId},
                        {fileSize: result.ContentLength,
                                location: data.Location,
                                status: "UPLOAD_COMPLETE"}).exec())
                    .then(async ()=>{
                        const url = `${process.env.CLOUD_RUN_URL}/run_classifier?uploadId=${uploadId}`;
                        const auth = new GoogleAuth();
                        const client = await auth.getIdTokenClient(url);
                        await client.request({url});
                        await projectModel.updateOne({_id: project!._id }, <any>{ status: ProjectJobStatus.PROCESSING_PENDING }, );

                        res.status(200).json({ project: project});
                    })
                    .catch((err)=>res.status(500).send(`Error persisting Multipart-Upload object data: ${err}`));
            });
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    })
    return router;
}