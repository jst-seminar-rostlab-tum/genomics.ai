import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";

export default function check_auth(){
    let router = express.Router();

    router.use((req : ExtRequest, res, next) => {
        if(req.header("auth") && req.header("auth") == "bearer"){
            // TODO LOGIC
            req.user_id = 42;
            req.isAuthenticated = true;
        }
    })

    return router;
}