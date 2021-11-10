import express, { Router } from "express";
import { tokenModel } from "../../../database/models/token";
import {userModel} from "../../../database/models/user";

export default function verify_email(): Router {
    let router = express.Router();

    router
        .get("/verify/:token", (async (req, res) => {
            const tokenObj =  await tokenModel.findOne({ token: req.params.token });

            if (!tokenObj)
                return res.status(404).send("Verification token could not be found. It may have expired.");
                
            const user = await userModel.findOne({ _id: tokenObj._userId });

            if (!user)
                return res.status(404).send("User for this token could not be found.");

            if (user.isEmailVerified) {
                res.status(200).send("User has already been verified.");
                return tokenObj.delete();
            }

            user.isEmailVerified = true;
            await user.save();
            res.status(200).send("User Email has been verified");
            tokenObj.delete();
        }))

    return router;
}