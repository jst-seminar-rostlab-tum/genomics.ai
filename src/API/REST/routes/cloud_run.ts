import express, {Router} from "express";
import { GoogleAuth } from "google-auth-library";
import { URL } from "url";

// Tests the Cloud Run connection
export default function cloud_run(): Router {
    let router = express.Router();

    router
        .post("/cloud_run", async (req: any, res) => {
            const {bucket_id} = req.body;

            if (!bucket_id)
                return res.status(400).send("Please provide a bucket ID!");

            const url = `${process.env.CLOUD_RUN_URL}?bucket_id=${bucket_id}`;

            const auth = new GoogleAuth();

            const client = await auth.getIdTokenClient(url);
            const response = await client.request({url});

            res.status(200).send(response.data);
        });

    return router;
}