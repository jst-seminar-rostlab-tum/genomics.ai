import express, { Router } from "express";
import { visibilityStatus } from "../../../../database/models/team";
import TeamService from "../../../../database/services/team.service";
import check_auth from "../../middleware/check_auth";
import { ExtRequest } from "../../../../definitions/ext_request";
import UserService from "../../../../database/services/user.service";

/**
 *  Returns all the teams that the user belong to.
 */
const get_teams_of_user = (): Router => {
  let router = express.Router();

  router.get("/users/ownteams", check_auth(), async (req: ExtRequest, res: any) => {
    try {
      const teams = req.user_id === undefined ? [] : await TeamService.getTeamsOfUser(req.user_id);
      return res.status(200).json(teams);
    } catch (e) {
      console.error("Error in get_teams_of_user()");
      console.error(JSON.stringify(e));
      console.error(e);
      return res.status(500).send("Internal error.");
    }
  });
  return router;
};

const get_users = (): Router => {
  let router = express.Router();

  router.get("/users", check_auth(), async (req: ExtRequest, res: any) => {
    try {
      const keyword = req.query?.keyword?.toString();
      const sort = req.query?.sort?.toString();
      const users = await UserService.searchUsers(keyword, sort);
      return res.status(200).json(users);
    } catch (e) {
      console.error("Error in /users");
      console.error(JSON.stringify(e));
      console.error(e);
      return res.status(500).send("Internal error.");
    }
  });
  return router;
};

const get_user_by_id = (): Router => {
  let router = express.Router();

  router.get("/users/:id", check_auth(), async (req: ExtRequest, res: any) => {
    try {
      const user_id = req.params.id;
      const {email, password, isEmailVerified, isAdministrator,
        ...filtered_user} = await UserService.getUserById(user_id) || {};

      if (!filtered_user) {
        return res.status(404).send("User with user id: " + user_id + " was not found!");
      }

      return res.status(200).json(filtered_user);
    } catch (e) {
      console.error("Error in /users/:id");
      console.error(JSON.stringify(e));
      console.error(e);
      return res.status(500).send("Internal error.");
    }
  });
  return router;
};

export { get_teams_of_user, get_users, get_user_by_id };
