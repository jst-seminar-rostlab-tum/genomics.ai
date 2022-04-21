import express, {Router} from "express";
import ProjectJobService from "../../../database/services/projectJob.service";
import { UpdateProjectJobDTO } from "../../../database/dtos/projectJob.dto";
import { GoogleAuth } from "google-auth-library";
import { ProjectJobStatus } from "../../../database/models/projectJob";
import check_auth from "../middleware/check_auth";

// Tests the Cloud Run connection
export default function abort_processing_route(): Router {
    let router = express.Router();

    router
        .post("/abort_processing",
            check_auth(),
            async (req: any, res) => {
                const {uploadId} = req.body;

                let project = await ProjectJobService.getProjectJobByUploadId(uploadId);

                if (!project)
                    return res.status(404).json({ msg: "No project found with upload ID." });

                if (project.owner != req.user_id)
                    return res.status(403).json({ msg: "A user may only abort their own projects!" });

                if (project.status != ProjectJobStatus[ProjectJobStatus.PROCESSING_PENDING])
                    return res.status(400).json({ msg: "Project processing cannot be aborted as it is not pending."});

                const update_object: UpdateProjectJobDTO = {
                    status: ProjectJobStatus[ProjectJobStatus.ABORTED]
                };
                await ProjectJobService.updateProjectJob(project._id, update_object);

                project.status = ProjectJobStatus[ProjectJobStatus.ABORTED];

                res.status(200).json({ project: project});
        });

    return router;
}
