import express, {Router} from "express";
import {userModel} from "../../../database/models/user";
import { tokenModel } from "../../../database/models/token";
import { send_mail } from "../../../util/mailer";

export default function approve_user(): Router {
    let router = express.Router();

    // This should be protected in the future so that only admins can invoke this.
    router
        .post("/approve", (async (req, res) => {
            /*const { academicAffiliation, email } = req.body;

            if (!email)
                return res.status(400).send('Please provide an e-mail address.');

            const user = await userModel.findOne({email});

            if (!user)
                return res.status(404).send('Found no user associated with the given e-mail address.');

            if (user.isAcademic)
                return res.status(200).send('User has already been approved.');

            user.isAcademic = true;
            user.academicAffiliation = academicAffiliation;
            await user.save();
            const token = await tokenModel.create({ _userId: user._id });

            send_mail(
                email, 
                'Your account has been approved',
                `Your account has been reviewed and approved by an administrator! Please verify here: http://localhost:8050/verify/${token.token}!`
            );

            res.status(200).send('User has been approved and verification e-mail has been sent.');
            */
        }));

    return router;
}