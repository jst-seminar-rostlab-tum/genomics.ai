import express from "express";
import {ExtRequest} from "../../../definitions/ext_request";
import {userModel} from "../../../database/models/user";
import check_auth from "../middleware/check_auth";

export default function get_profile_route() {
    let router = express.Router();

    router.get("/profile", check_auth(), async (req: ExtRequest, res: any) => {
        const user = await userModel.findById(req.user_id);

        return res.status(200).json(user);
    });
    return router;
}
