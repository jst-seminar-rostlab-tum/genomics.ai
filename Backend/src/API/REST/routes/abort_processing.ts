import express, { Router } from "express";
import ProjectService from "../../../database/services/project.service";
import { UpdateProjectDTO } from "../../../database/dtos/project.dto";
import { GoogleAuth } from "google-auth-library";
import { ProjectStatus } from "../../../database/models/project";
import check_auth from "../middleware/check_auth";
import { validationMdw } from "../middleware/validation";

// Tests the Cloud Run connection
export default function abort_processing_route(): Router {
  let router = express.Router();

  router.post("/abort_processing", validationMdw, check_auth(), async (req: any, res) => {
    const { uploadId } = req.body;

    let project = await ProjectService.getProjectByUploadId(uploadId);

    if (!project) return res.status(404).json({ msg: "No team found with upload ID." });

    if (project.owner != req.user_id)
      return res.status(403).json({ msg: "A user may only abort their own projects!" });

    if (project.status != ProjectStatus[ProjectStatus.PROCESSING_PENDING])
      return res
        .status(400)
        .json({ msg: "Project processing cannot be aborted as it is not pending." });

    const update_object: UpdateProjectDTO = {
      status: ProjectStatus[ProjectStatus.ABORTED],
    };
    await ProjectService.updateProjectById(project._id, update_object);

    project.status = ProjectStatus[ProjectStatus.ABORTED];

    res.status(200).json({ project: project });
  });

  return router;
}
