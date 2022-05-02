import express, { Router } from "express";
import UserService from "../../../database/services/user.service";
import PasswordResetTokenService from "../../../database/services/passwordResetToken.service";
import { AddPasswordResetTokenDTO } from "../../../database/dtos/password_reset_token.dto";
import bcrypt from "bcrypt";
import { mailer } from "../../../util/mailer";
import { validationMdw } from "../middleware/validation";

export default function password_reset_route(): Router {
  let router = express.Router();

  router.post("/password_reset", validationMdw, async (req, res) => {
    const { email } = req.body;
    if (!email) return res.status(400).send("Missing parameters");

    // check user with this email exists
    const user = await UserService.getUserByEmail(email);
    if (!user) return res.status(404).send("User not found");

    // TODO: check for existing tokens - don't send mail if there already is a recent token. Delete old tokens!

    // create a reset-token and send email
    const tokenToAdd: AddPasswordResetTokenDTO = {
      _userId: user._id,
    };
    const token = await PasswordResetTokenService.addToken(tokenToAdd);
    await mailer.send(
      user.email,
      "[GeneCruncher] Please reset your password",
      "password_reset_request_email",
      {
        firstname: user.firstName,
        // link: `https://www.genecruncher.com/#/password_reset?token=${token.token}`
        link: `${process.env.FRONTEND_URL}/#/password_reset?token=${token.token}`,
        new_reset_link: `${process.env.FRONTEND_URL}/password_reset`,
      }
    );

    return res
      .status(200)
      .send(
        "An email has been sent to your email address with instructions to reset your password"
      );
  });

  router.post("/password_reset/:token", validationMdw, async (req, res) => {
    const { password } = req.body;
    if (!password) return res.status(400).send("Missing parameters");

    // find token corresponding to token-id and delete it
    const token = await PasswordResetTokenService.getTokenByToken(req.params.token);
    if (!token)
      return res
        .status(404)
        .send("Token could not be found. It may have been expired or used already");
    token.deleteOne();

    // find user indicated by token's user-field
    const user = await UserService.getUserById(token._userId);
    if (!user)
      return res.status(404).send("Invalid token: User for this token could not be found.");

    // update the password
    user.password = await bcrypt.hash(password, 15);
    user.save();

    // send pw-update email
    await mailer.send(
      user.email,
      "[GeneCruncher] Your password was successfully reset",
      "password_reset_confirmation_email",
      {
        firstname: user.firstName,
        email: user.email,
        link: `${process.env.FRONTEND_URL}/password_reset`,
      }
    );

    return res.status(200).send("Password has been changed");
  });

  return router;
}
