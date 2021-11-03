import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import {projectModel} from "../../../database/models/project";
import check_auth from "../middleware/check_auth";

export default function get_job_route() {
  let router = express.Router();

  router
    .get("/get_jobs", check_auth(), async (req: ExtRequest, res: any) => {
      try {
        const jobs = await projectModel.find({owner: req.user_id});
        return res.status(200).json(jobs!);
      } catch (e) {
        return res.status(500);
      }
    });
  return router;
}
