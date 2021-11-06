import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import {userModel} from "../../../database/models/user";
import check_auth from "../middleware/check_auth";

export default function get_profile_route() {
  let router = express.Router();

  router
    .get("/profile/:id", check_auth(), async (req: ExtRequest, res: any) => {
      const userId = req.params.id;

      const user = await userModel.findById(userId);

      if (!user) {
        res.status(404).send(`User ${userId} not found`);
      }

      return res.status(200).json(user);
    });
  return router;
}
