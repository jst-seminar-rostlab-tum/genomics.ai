import express, {Router} from "express";
import { GoogleAuth } from "google-auth-library";
import { IJob, jobModel } from "../../../database/models/job";

// Tests the Cloud Run connection
export default function cloud_run(): Router {
    let router = express.Router();

    router
        .post("/cloud_run", async (req: any, res) => {
            const {bucket_id} = req.body;

            if (!bucket_id)
                return res.status(400).send("Please provide a bucket ID!");

            if (await jobModel.findOne({bucketId: bucket_id}))
                return res.status(400).send("There already exists a job for the given bucket ID!");

            let job: (IJob | undefined) = undefined;
            try {
                job = await jobModel.create({
                    bucketId: bucket_id
                });
            }
            catch(err) {
                return res.status(500).json({ msg: "Unable to create job."});
            }

            const url = `${process.env.CLOUD_RUN_URL}?bucket_id=${bucket_id}`;

            const auth = new GoogleAuth();

            const client = await auth.getIdTokenClient(url);
            const response = await client.request({url});

            res.status(200).json({ msg: response.data, job: job});
        });

    return router;
}