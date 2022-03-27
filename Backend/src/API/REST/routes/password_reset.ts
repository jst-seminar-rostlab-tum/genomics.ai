import express, {Router} from "express";
import {userModel} from "../../../database/models/user";
import {passwordResetTokenModel} from "../../../database/models/password_reset_token";
import bcrypt from "bcrypt";
import {mailer} from "../../../util/mailer";

export default function password_reset_route() : Router {
    let router = express.Router();
    
    router.post('/password_reset', async (req, res) => {
      const {email} = req.body;
      if (!email)
        return res.status(400).send('Missing parameters');

      // check user with this email exists
      const user = await userModel.findOne({email});
      if (!user)
        return res.status(404).send("User not found");

      // TODO: check for existing tokens - don't send mail if there already is a recent token. Delete old tokens!

      // create a reset-token and send email
      const token = await passwordResetTokenModel.create({
        _userId: user._id,
      });
      await mailer.send(user.email, "[GeneCruncher] Please reset your password", "password_reset_request_email", {
        firstname: user.firstName,
        link: `https://www.genecruncher.com/#/password_reset?token=${token.token}`
      })

      return res.status(200).send('An email has been sent to your email address with instructions to reset your password');
    });
    
    router.post('/password_reset/:token', async (req, res) => {
      const {password} = req.body;
      if (!password)
        return res.status(400).send("Missing parameters");

      // find token corresponding to token-id and delete it
      const token = await passwordResetTokenModel.findOne({token: req.params.token});
      if (!token)
        return res.status(404).send("Token could not be found. It may have been expired or used already");
      token.deleteOne();

      // find user indicated by token's user-field
      const user = await userModel.findById(token._userId);
      if (!user)
        return res.status(404).send("Invalid token: User for this token could not be found.");

      // update the password
      user.password = await bcrypt.hash(password, 15);
      user.save();

      // send pw-update email
      await mailer.send(user.email, "[GeneCruncher] Your password was successfully reset", "password_reset_confirmation_email", {
        firstname: user.firstName,
        email: user.email
      });

      return res.status(200).send("Password has been changed");
    });
    
    return router;
}
