import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import jwt from "jsonwebtoken";
import {userModel} from "../../../database/models/user";

export default function check_auth(){
    let router = express.Router();

    router.use((req : ExtRequest, res, next) => {
        req.is_authenticated = false;
        if(req.header("auth")){
            // TODO: Secret => REDIS
            try{
                jwt.verify(req.header("auth")!, "SECRET", async function(err, decoded){
                    if(err || !decoded || !decoded.email){
                        if(err?.name == "TokenExpiredError")
                            return res.status(440).send("JWT authentication token expired. Please log in again")

                        return res.status(401).send("Invalid authentication");
                    }

                    userModel.findOne({_id: decoded.id}).exec((err, result)=>{
                        if(err){
                            console.error(err);
                            return res.status(500).send("Error during authentication: Failed to fetch user")
                        }
                        req.is_authenticated = true;
                        req.user_id = decoded.id;
                        req.email = result!.email;
                        next();
                    })
                });
            }catch(e){
                return console.error(e); // abort on error;
            }
        }
    })

    return router;
}