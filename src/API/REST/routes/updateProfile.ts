import express from "express";
import verify_authenticated_user from "../middleware/verify_authenticated_user";
import {ExtRequest} from "../../../definitions/ext_request";

export default function update_profile_route(){
    let router = express.Router();

    router
        .post(  "/update_profile",
                verify_authenticated_user(),
                (req : ExtRequest, res) => {
            console.log(`${req.user_id} has updated their profile`);
            res.send("Your profile has been updated")
        })

    return router;
}