import express from "express";
import { userModel } from "../../../database/models/user";
import { ExtRequest } from "../../../definitions/ext_request";

export default function verify_email_route() {
  let router = express.Router();

  router
    .post("/verify_email/:token", async (req: ExtRequest, res: any) => {
      const verificationToken = req.params.token;

      try {
        await userModel.findOneAndUpdate({emailVerificationToken: verificationToken}, {isVerified: true});

        return res.status(200).send("Your email address has been verified!");
      } catch (err) {
        return res.status(500).send(err);
      }
    });
  return router;
}
