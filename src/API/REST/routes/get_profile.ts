import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import {userModel} from "../../../database/models/user";
import check_auth from "../middleware/check_auth";

export default function get_profile_route() {
  let router = express.Router();

  router
    .get("/get_profile/:id", check_auth(), async (req: ExtRequest, res: any) => {
      const userId = req.params.id;

      try {
        const user = await userModel.findOne({_id: userId});
        return res.status(200).json(user!.toJSON());
      } catch (e) {
        return res.status(500);
      }
    });
  return router;
}
