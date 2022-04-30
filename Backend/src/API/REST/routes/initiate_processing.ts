import express, { Router } from "express";
import { GoogleAuth } from "google-auth-library";
import { ProjectStatus } from "../../../database/models/project";
import ProjectService from "../../../database/services/project.service";
import check_auth from "../middleware/check_auth";

// Tests the Cloud Run connection
export default function initiate_processing_route(): Router {
  let router = express.Router();
  //DISABLED AT THE MOMENT probably never used? Seems to be part of the test routes?
  router.post("/initiate_processing", check_auth(), async (req: any, res) => {
    const { uploadId } = req.body;

    let project = await ProjectService.getProjectByUploadId(uploadId);

    if (!project) return res.status(404).json({ msg: "No team found with upload ID." });

    if (project.owner != req.user_id)
      return res.status(403).json({ msg: "A user may only initiate their own projects!" });

    if (project.status != ProjectStatus[ProjectStatus.UPLOAD_COMPLETE])
      return res.status(400).json({
        msg: "Processing cannot be initiated. The upload has to be finished uploading and can only be initiated once.",
      });

    const url = `${process.env.CLOUD_RUN_URL}/run_classifier?uploadId=${uploadId}`;

    const auth = new GoogleAuth();

    const client = await auth.getIdTokenClient(url);
    await client.request({ url });

    res.status(200).json({ project: project });
  });

  return router;
}
