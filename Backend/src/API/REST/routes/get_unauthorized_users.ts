import express from "express";
import UserService from "../../../database/services/user.service";
import { ExtRequest } from "../../../definitions/ext_request";
import check_auth from "../middleware/check_auth";

export default function get_unauthorized_users_route() {
  let router = express.Router();

  router.get("/unauthorized_users", check_auth(), async (req: ExtRequest, res) => {
    if (!req.is_administrator) return res.status(403).send("Unauthorized");

    try {
      const unauthorizedUsers = await UserService.getUnauthUsers();
      return res.status(200).json(unauthorizedUsers!);
    } catch (err) {
      return res.status(500).send(err);
    }
  });
  return router;
}
