import express from "express";
import jwt from "jsonwebtoken";
import {userModel} from "../../../database/models/user";
import bcrypt from "bcrypt";

const INCORRECT_CREDENTIALS =  "The email or password is incorrect"; 
const JWT_SECRET = process.env.JWT_SECRET || ""; 

export default function auth_route(){
    let router = express.Router();

    router.post("/auth", (req, res, next) => {
        const {email, password} = req.body;
        userModel.findOne({email: email}).select('+password').exec().then(user =>{
            if(!user)
                return res.status(401).send(INCORRECT_CREDENTIALS);
            if(!user.isEmailVerified)
                return res.status(401).send("User not verified");

            bcrypt.compare(password, <string>user.password, (err, match) => {
                if (err) {
                    console.error(err);
                    return res.status(500).send(INCORRECT_CREDENTIALS);
                }

                if (match) {
                    delete user.password;
                    const token = jwt.sign(
                        {id: user._id, email: user.email},
                        JWT_SECRET,
                        { expiresIn: '20h' });

                    /* user without the password field */
                    const { password, ...userSecure } = user.toObject();

                    if(!userSecure.avatarUrl) {
                        userSecure.avatarUrl=`${process.env.S3_ENDPOINT}/${process.env.S3_PICTURE_BUCKET_NAME}/user.png`
                    }

                    return res.status(200).json({
                        msg: "Login success",
                        user: userSecure,
                        jwt: token
                    })
                } else {
                    return res.status(401).json({ msg: INCORRECT_CREDENTIALS })
                }
            })
        })
    })

    return router;
}
