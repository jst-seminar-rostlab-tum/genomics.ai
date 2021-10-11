import express, { Router } from "express";
import { tokenModel } from "../../../database/models/token";
import {userModel} from "../../../database/models/user";
import { send_verification_mail } from "../../../util/mailer";

export default function resend_verification_link(): Router {
    let router = express.Router();

    router
        .post("/resend", (async (req, res) => {
            const email = req.body.email;

            if (!email)
                return res.status(400).send("Missing e-mail parameter.");

            const user = await userModel.findOne({email});

            if (!user)
                return res.status(404).send("Could not find user with this e-mail. Please register.");

            if (user.isVerified)
                return res.status(200).send("User has already been verified.");

            let token = await tokenModel.findOne({_userId: user._id});

            if (token)
                await token.delete();

            token = await tokenModel.create({ _userId: user._id });

            send_verification_mail(email, token.token);
            res.status(200).send("Verification link has been resent.");
        }));

    return router;
}