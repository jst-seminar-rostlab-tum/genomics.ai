import express from "express";
import { ExtRequest } from "../../../definitions/ext_request";
import UserService from "../../../database/services/user.service";
import { UpdateUserDTO } from "../../../database/dtos/user.dto";
import bcrypt from "bcrypt";
import check_auth from "../middleware/check_auth";
import { validationMdw } from "../middleware/validation";

export default function update_profile_route() {
  let router = express.Router();

  router.post("/update_profile", validationMdw, check_auth(), async (req: ExtRequest, res) => {
    const { first_name, last_name, email, password, note } = req.body;

    if (!(first_name || last_name || email || password || note))
      return res.status(400).send("Missing parameters");

    let userFound = req.user_id === undefined ? false : await UserService.getUserById(req.user_id);
    if (!userFound) return res.status(404).send("User not found");

    try {
      let update_object: UpdateUserDTO = {};
      // TODO check if they are different
      if (first_name) update_object.firstName = first_name;
      if (last_name) update_object.lastName = last_name;
      if (note) update_object.note = note;
      if (email)
        // TODO implement changing email addresses
        return res.status(501).send("Changing email-address is not implemented yet.");
      if (password) update_object.password = await bcrypt.hash(password, 15);

      if (req.user_id !== undefined) await UserService.updateUser(req.user_id, update_object);

      res.status(200).send({ msg: "User updated." });
    } catch (err) {
      console.error(err);
      return res.status(500).send("Unable to create user. (DB-error)");
    }
  });

  return router;
}
