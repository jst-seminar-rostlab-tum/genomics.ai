import express, { Router } from "express";
import bcrypt from "bcrypt";
import { IUser, userModel } from "../../../database/models/user";

export default function hello_route(): Router {
  let router = express.Router();

  router.post("/hello", async (req: any, res) => {
    // TODO input-validation
    const { ping } = req.body;

    let rand_float = Math.random();

    if (rand_float > 0.75) return res.status(500).send("Some server error occurred");
    if (rand_float > 0.5) return res.status(400).send("Some user error occurred");

    return res.status(200).json({ msg: "Pong: " + ping });
  });

  return router;
}
