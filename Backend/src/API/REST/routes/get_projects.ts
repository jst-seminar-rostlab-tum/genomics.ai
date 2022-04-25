import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import ProjectService from "../../../database/services/project.service";
import check_auth from "../middleware/check_auth";

export default function get_projects_route() {
    let router = express.Router();

    router
    .get("/projects", check_auth(), async (req: ExtRequest, res: any) => {
        try {
            const projects = req.user_id === undefined ? null :
              await ProjectService.getProjectByOwner(req.user_id, -1);

            return res.status(200).json(projects!);
        } catch (e) {
            return res.status(500);
        }
    });
    return router;
}
