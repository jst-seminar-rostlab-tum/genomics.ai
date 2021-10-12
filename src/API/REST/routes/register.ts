import express, {Router} from "express";
import {userModel} from "../../../database/models/user";
import { tokenModel } from "../../../database/models/token";
import bcrypt from "bcrypt";
import { send_verification_mail } from "../../../util/mailer";
import { Verifier } from 'academic-email-verifier';

export default function register_route(): Router {
    let router = express.Router();

    router
        .post("/register", (async (req, res) => {
            const {firstName, lastName, email, password} = req.body;

            if (await userModel.findOne({email})) {
                return res.status(409).send("User with the given email already exists");
            }

            const encryptedPassword = await bcrypt.hash(password, 15);
            const institutionName = await Verifier.getInstitutionName(email);
            const isAcademic = institutionName !== "";
            
            const user = await userModel.create({
                firstName, lastName, email, 
                password: encryptedPassword,
                academicAffiliation: institutionName,
                isAcademic: isAcademic
            });

            if (isAcademic) {
                const token = await tokenModel.create({ _userId: user._id });
                send_verification_mail(email, token.token);
                return res.status(201).json(user);
            }

            return res.status(403).send("The e-mail address does not seem to belong to an academic institution. Wait for an administrator to manually approve you.");
        }))

    return router;
}