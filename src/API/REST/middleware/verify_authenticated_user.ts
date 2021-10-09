import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";

export default function verify_authenticated_user(){
    let router = express.Router();

    router.use((req : ExtRequest, res, next) => {
        if(req.header("auth") && req.header("auth") == "bearer"){
            // TODO LOGIC
            req.isAuthenticated = true;
        }else{
            return res
                .status(401)
                .send("Missing authentication");
        }

        next();
    })

    return router;
}