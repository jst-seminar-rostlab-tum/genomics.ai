import express, { Router } from "express";
import bcrypt from "bcrypt";
import { IUser, userModel } from "../../../database/models/user";
import check_auth from "../middleware/check_auth";

export default function hello_auth_route(): Router {
  let router = express.Router();

  router.post("/hello", check_auth(), async (req: any, res) => {
    // TODO input-validation
    const { ping } = req.body;

    userModel
      .findOne({ _id: req.user_id })
      .exec()
      .then((user) => {
        if (!user) return res.status(401).send("User not found");

        return res.status(200).json({ msg: `Hello ${user.firstName}! Pong: ${ping}` });
      });
  });

  return router;
}
