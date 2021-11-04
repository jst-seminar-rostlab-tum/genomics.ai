import express, {Router} from "express";
import bcrypt from "bcrypt";
import crypto from "crypto";
import {IUser, userModel} from "../../../database/models/user";

export default function register_route(): Router {
    let router = express.Router();

    router
        .post("/register", async (req: any, res) => {
            // TODO input-validation
            const {first_name, last_name, email, password, note} = req.body;

            if (!(first_name && email && password))
                return res.status(400).send("Missing parameters");

            if (await userModel.findOne({email}))
                return res.status(409).send("User with the given email already exists");

            const saltHashedPassword = await bcrypt.hash(password, 15);
            const emailVerificationToken = await crypto.randomBytes(16).toString("hex");

            let user : (IUser | undefined) = undefined;
            try{
                await userModel.create({
                    firstName: first_name, lastName: last_name, email, password: saltHashedPassword, note, emailVerificationToken
                });

                return res
                    .status(201)
                    .json({
                        msg: "User has been registered"
                    })
            }catch(err){
                console.error(err);
                return res.status(500).send("Unable to create user. (DB-error)");
            }
        })

    return router;
}
