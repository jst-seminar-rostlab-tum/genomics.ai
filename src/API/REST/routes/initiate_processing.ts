import express, {Router} from "express";
import { GoogleAuth } from "google-auth-library";
import { projectModel, ProjectJobStatus } from "../../../database/models/project";
import check_auth from "../middleware/check_auth";
import { mongo } from "mongoose";

// Tests the Cloud Run connection
export default function initiate_processing_route(): Router {
    let router = express.Router();

    router
        .post("/initiate_processing",
            check_auth(),
            async (req: any, res) => {
                const {uploadId} = req.body;

                let project = await projectModel.findOne({ uploadId: uploadId });

                if (!project)
                    return res.status(404).json({ msg: "No project found with upload ID." });

                if (project.owner != req.user_id)
                    return res.status(403).json({ msg: "A user may only initiate their own projects!" });

                if (project.status != ProjectJobStatus.UPLOAD_COMPLETE.toString())
                    return res.status(400).json({ msg: "Processing cannot be initiated. The upload has to be finished uploading and can only be initiated once."});

                const url = `${process.env.CLOUD_RUN_URL}?upload_id=${uploadId}`;

                const auth = new GoogleAuth();

                const client = await auth.getIdTokenClient(url);
                await client.request({url});

                await projectModel.updateOne({_id: project._id }, <any>{ status: ProjectJobStatus.PROCESSING_PENDING }, );
                project.status = ProjectJobStatus.PROCESSING_PENDING.toString();
                
                res.status(200).json({ project: project});
        });

    return router;
}