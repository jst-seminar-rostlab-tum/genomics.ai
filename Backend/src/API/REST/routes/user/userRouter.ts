import express, {Router} from "express";
import { visibilityStatus } from "../../../../database/models/team";
import TeamService from "../../../../database/services/team.service";
import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";

/**
 *  Returns all the teams of the user.
 */
const get_teams_of_user = (): Router => {
    let router = express.Router();

    router.get("/users/ownteams", check_auth(), async (req: ExtRequest, res: any) => {
        try {
            const teams = req.user_id === undefined ? [] :
              await TeamService.getTeamsOfUser(req.user_id);
            return res.status(200).json(teams);
        } catch(e) {
            console.error("Error in get_teams_of_user()");
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }
    })
    return router;
}

export { get_teams_of_user };
