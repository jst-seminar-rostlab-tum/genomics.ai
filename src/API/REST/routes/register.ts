import express, {Router} from "express";
import bcrypt from "bcrypt";
import {IUser, userModel} from "../../../database/models/user";
import {Verifier} from "academic-email-verifier";
import {send_verification_mail} from "../../../util/mailer";
import {tokenModel} from "../../../database/models/token";

export default function register_route() : Router {
    let router = express.Router();

    router
        .post("/register", async (req: any, res) => {
            // TODO input-validation
            const {first_name, last_name, email, password, note} = req.body;

            if (!(first_name && email && password))
                return res.status(400).send("Missing parameters");

            if (await userModel.findOne({email}))
                return res.status(409).send("User with the given email already exists");

            const institutionName = await Verifier.getInstitutionName(email);
            const isAcademic = institutionName !== "";
            const saltHashedPassword = await bcrypt.hash(password, 15);

            let user : (IUser | undefined) = undefined;
            try{
                user = await userModel.create({
                    firstName: first_name,
                    lastName: last_name,
                    email,
                    password: saltHashedPassword,
                    note,
                    isAccountApproved: isAcademic
                });

                if (isAcademic) {
                    const token = await tokenModel.create({ _userId: user._id });
                    send_verification_mail(email, token.token);
                    return res.status(201).json(user);
                }else{
                    return res.status(200).send("The e-mail address does not seem to belong to an academic institution. Wait for an administrator to manually approve you.");
                }

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