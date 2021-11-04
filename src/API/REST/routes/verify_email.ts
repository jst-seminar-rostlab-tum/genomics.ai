import express from "express";
import { userModel } from "../../../database/models/user";
import { ExtRequest } from "../../../definitions/ext_request";

export default function verify_email_route() {
  let router = express.Router();

  router
    .get("/verify_email/:token", async (req: ExtRequest, res: any) => {
      const verificationToken = req.params.token;

      try {
        await userModel.findOneAndUpdate({emailVerificationToken: verificationToken}, {isVerified: true});

        return res
          .set("Content-Type", "text/html")
          .send(Buffer.from("<h1>Email verified</h2><p>Congratulations, you have successfully verified your email address!</p>"));
      } catch (err) {
        return res.status(500).send(err);
      }
    });
  return router;
}
