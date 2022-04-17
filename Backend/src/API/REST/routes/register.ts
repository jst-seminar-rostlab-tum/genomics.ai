import express, {Router} from "express";
import bcrypt from "bcrypt";
import {IUser} from "../../../database/models/user";
import {AddUserDTO} from "../../../database/dtos/user.dto"
import UserService from "../../../database/services/user.service";
import {mailer} from "../../../util/mailer";
import {tokenModel} from "../../../database/models/token";

export default function register_route() : Router {
    let router = express.Router();

    router
        .post("/register", async (req: any, res) => {
            const {first_name, last_name, email, password, note} = req.body;

            if (!(first_name && email && password))
                return res.status(400).send("Missing parameters");

            if (await UserService.getUserByEmail(email))
                return res.status(409).send("User with the given email already exists");

            const saltHashedPassword = await bcrypt.hash(password, 12);

            let userToAdd: AddUserDTO = {
                firstName: first_name,
                lastName: last_name,
                email,
                password: saltHashedPassword,
                note,
            }

            let userAdded: (IUser | undefined) = undefined;

            try {
                userAdded = await UserService.addUser(userToAdd);

                /* user without the password field */
                const { password, ...userSecure } = userAdded.toObject();

                const token = await tokenModel.create({ _userId: userAdded._id });
                mailer.send_verification_mail(first_name, email, token.token);
                return res.status(201).json(userSecure);
            } catch(err) {
                console.error("Error registering user!");
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to create user. (DB-error)");
            }
        })

    return router;
}
