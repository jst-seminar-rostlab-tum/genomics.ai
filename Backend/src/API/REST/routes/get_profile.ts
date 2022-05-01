import express from "express";
import { ExtRequest } from "../../../definitions/ext_request";
import UserService from "../../../database/services/user.service";
import check_auth from "../middleware/check_auth";

export default function get_profile_route() {
  let router = express.Router();

  router.get("/profile", check_auth(), async (req: ExtRequest, res: any) => {
    const user = req.user_id === undefined ? null : await UserService.getUserById(req.user_id);

    return res.status(200).json(user);
  });
  return router;
}
