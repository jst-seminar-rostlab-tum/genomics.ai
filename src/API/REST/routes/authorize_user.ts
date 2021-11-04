import express from "express";
import { userModel } from "../../../database/models/user";
import { ExtRequest } from "../../../definitions/ext_request";
import check_auth from "../middleware/check_auth";

export default function authorize_user_route() {
  let router = express.Router();

  router
    .post("/authorize_user/:id", check_auth(), async (req: ExtRequest, res: any) => {
      if (!req.is_administrator) {
        return res.status(403).send("Unauthorized");
      }

      const userId = req.params.id;

      try {
        await userModel.findByIdAndUpdate(userId, {authorized: true});
        return res.status(200).send("User has been authorized");
      } catch (err) {
        return res.status(500).send(err);
      }
    });
  return router;
}
