import express, {Router} from "express";
import { GoogleAuth } from "google-auth-library";
import { URL } from "url";

// Tests the Cloud Run connection
export default function cloud_run(): Router {
    let router = express.Router();

    router
        .post("/cloud_run", async (req: any, res) => {
            const {bucket_id} = req.body;
            const url = <string>process.env.CLOUD_RUN_URL;

            const auth = new GoogleAuth();
            const targetAudience = new URL(url);

            const client = await auth.getIdTokenClient(url);
            const response = await client.request({url});

            res.status(200).send(response.data);
        });

    return router;
}