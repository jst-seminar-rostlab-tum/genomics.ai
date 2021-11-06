import express from "express";
import { userModel } from "../../../database/models/user";
import { ExtRequest } from "../../../definitions/ext_request";

export default function verify_email_route() {
  let router = express.Router();

  router
    .get("/verify_email/:id", async (req: ExtRequest, res: any) => {
      const userId = req.params.id;
      const verificationToken = req.query.token;

      if (verificationToken == null) {
        return res.status(400).send("Missing token query parameter");
      }

      const user = await userModel.findById(userId);

      if (!user) {
        return res.status(404).send(`User ${userId} not found`);
      }

      if (user.emailVerificationToken !== verificationToken) {
        return res.status(400).send("Wrong token");
      }

      user.isVerified = true;
      await user.save();

      return res
        .set("Content-Type", "text/html")
        .send(Buffer.from("<h1>Email verified</h2><p>Congratulations, you have successfully verified your email address!</p>"));
    });
  return router;
}
