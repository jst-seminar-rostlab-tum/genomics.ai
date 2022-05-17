import express from "express";
import UserService from "../../../database/services/user.service";
import { ExtRequest } from "../../../definitions/ext_request";
import check_auth from "../middleware/check_auth";

// TODO: This route requires JWT in cookies to be set (for the browser to send it with any GET-request automatically)

export default function authorize_user_route() {
  let router = express.Router();

  router.get("/authorize_user/:id", check_auth(), async (req: ExtRequest, res: any) => {
    if (!req.is_administrator) return res.status(403).send("Unauthorized");

    const userId = req.params.id;
    const user = await UserService.getUserById(userId);

    if (!user) return res.status(404).send(`User ${userId} not found`);

    // set authorized to true
    user.isAuthorized = true;
    user.isEmailVerified = true;
    await user.save();

    return res.status(200).send(`User ${userId} has been authorized`);
  });
  return router;
}
