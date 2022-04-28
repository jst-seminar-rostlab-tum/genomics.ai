import check_auth from "../../middleware/check_auth";
import { ExtRequest } from "../../../../definitions/ext_request";
import ProjectService from "../../../../database/services/project.service";
import { UploadPartRequest } from "aws-sdk/clients/s3";
import s3 from "../../../../util/s3";
import express from "express";

export default function upload_get_upload_url_route() {
    let router = express.Router();

    router.get("/file_upload/get_upload_url", check_auth(), async (req: ExtRequest, res) => {
        let { partNumber, uploadId } = req.query;
        if (!process.env.S3_BUCKET_NAME) return res.status(500).send("S3-BucketName is not set");

        try {
            let project = await ProjectService.getProjectByUploadId(String(uploadId), req.user_id);

            if (!project) return res.status(400).send("Upload not found");

            let params: UploadPartRequest = {
                Bucket: process.env.S3_BUCKET_NAME!,
                Key: String(project._id),
                PartNumber: Number(partNumber),
                UploadId: String(uploadId)
            };
            let presignedUrl = await s3.getSignedUrlPromise("uploadPart", params);
            res.status(200).send({ presignedUrl });
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    });
    return router;
}
