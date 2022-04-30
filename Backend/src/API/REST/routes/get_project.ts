import express from "express";
import { ExtRequest } from "../../../definitions/ext_request";
import ProjectService from "../../../database/services/project.service";
import check_auth from "../middleware/check_auth";

export default function get_project_route() {
  let router = express.Router();

  router.get("/project/:id", check_auth(), async (req: ExtRequest, res: any) => {
    const projectId = req.params.id;
    const project = await ProjectService.getProjectById(projectId);

    if (!project) return res.status(404).send(`Project ${projectId} not found`);

    return res.status(200).json(project);
  });
  return router;
}
