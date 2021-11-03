import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import {projectModel} from "../../../database/models/project";
import check_auth from "../middleware/check_auth";

export default function get_job_route() {
  let router = express.Router();

  router
    .get("/get_job/:id", check_auth(), async (req: ExtRequest, res: any) => {
      const jobId = req.params.id;

      try {
        const job = await projectModel.findOne({_id: jobId});
        return res.status(200).json(job!.toJSON());
      } catch (e) {
        return res.status(500);
      }
    });
  return router;
}
