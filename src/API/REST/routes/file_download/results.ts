import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";
import {ProjectJobStatus, projectModel} from "../../../../database/models/project";
import s3 from "../../../../util/s3";
import express from "express";

export default function download_results_route() {
    let router = express.Router();
    router.post('/file_download/results', check_auth(), async (req: ExtRequest, res) => {
        let {id} = req.body;
        if(!id)
            return res.status(400).send("Missing job id parameter.");

        try {
            if (process.env.S3_BUCKET_NAME && req.user_id) {
                const job = await projectModel.findById(id);
                if (!job) {
                    return res.status(404).send("Job not found.");
                }
                if (job.status != ProjectJobStatus[ProjectJobStatus.DONE]) {
                    return res.status(400).send("Job not completed.");
                }
                if (job.owner != req.user_id) {
                    return res.status(400).send("User not authorized to access file.");
                }

                let params: any = {
                    Bucket: process.env.S3_BUCKET_NAME!,
                    Key: String(job.resultName),
                }
                let presignedUrl = await s3.getSignedUrlPromise('getObject', params);
                res.status(200).send({presignedUrl});
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    })
    return router;
}