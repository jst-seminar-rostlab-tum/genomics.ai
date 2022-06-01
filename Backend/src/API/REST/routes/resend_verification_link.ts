import express, { Router } from "express";
import { AddTokenDTO } from "../../../database/dtos/token.dto";
import TokenService from "../../../database/services/token.service";
import UserService from "../../../database/services/user.service";
import { mailer } from "../../../util/mailer";
import { validationMdw } from "../middleware/validation";

export default function resend_verification_link(): Router {
  let router = express.Router();

  router.post("/resend", validationMdw, async (req, res) => {
    const email = req.body.email;
    if (!email) return res.status(400).send("Missing e-mail parameter.");

    // check if user needs verification
    const user = await UserService.getUserByEmail(email);
    if (!user)
      return res.status(404).send("Could not find user with this e-mail. Please register.");
    if (user.isEmailVerified) return res.status(200).send("User has already been verified.");

    // delete old token
    const token = await TokenService.getTokenByUserId(user._id);
    // TODO verify token is old enough (has to be at least 60 seconds old)
    if (token) await token.delete();

    // create new token and send new verification-mail
    const tokenToAdd: AddTokenDTO = { _userId: user._id };
    const tokenNew = await TokenService.addToken(tokenToAdd);
    await mailer.send_verification_mail(user.firstName, email, tokenNew.token);

    res.status(200).send("Verification link has been resent.");
  });

  return router;
}
