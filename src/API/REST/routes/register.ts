import express, {Router} from "express";
import bcrypt from "bcrypt";

import {userModel} from "../../../database/models/user";

export default function register_route(): Router {
    let router = express.Router();

    router
        .post("/register")
        .use((async (req, res) => {
            // TODO input validation
            const {firstName, lastName, email, password} = req.body;

            if (await userModel.findOne({email})) {
                return res.status(409).send("User with the given email already exists");
            }

            const encryptedPassword = await bcrypt.hash(password, 15);

            const user = await userModel.create({
                firstName, lastName, email, password: encryptedPassword
            })

            res.status(201).json(user);
        }))

    return router;
}