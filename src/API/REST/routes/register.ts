import express, {Router} from "express";
import bcrypt from "bcrypt";
import jwt from "jsonwebtoken";
import {userModel} from "../../../database/models/user";
import {ExtRegisterRequest} from "../../../definitions/ext_register_request";

export default function register_route(): Router {
    let router = express.Router();

    router
        .post("/register")
        .use((async (req: ExtRegisterRequest, res) => {
            const {firstName, lastName, email, password} = req.body;

            if (!(firstName && lastName && email && password)) {
                return res.status(400).send("Missing parameters");
            }

            if (await userModel.findOne({email})) {
                return res.status(409).send("User with the given email already exists");
            }

            const encryptedPassword = await bcrypt.hash(password, 15);

            const user = await userModel.create({
                firstName, lastName, email, password: encryptedPassword
            });

            // TODO: Secret
            const token = jwt.sign({id: user._id, email: user.email}, "SECRET");

            return res
                .cookie("access_token", token, {httpOnly: true, secure: true})
                .status(201)
                .send("User has been registered")
        }))

    return router;
}