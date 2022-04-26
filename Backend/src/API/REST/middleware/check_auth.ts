import express from "express";
import { ExtRequest } from "../../../definitions/ext_request";
import jwt from "jsonwebtoken";
import UserService from "../../../database/services/user.service";

const JWT_SECRET = process.env.JWT_SECRET || "";

export default function check_auth() {
  let router = express.Router();

  router.use((req: ExtRequest, res, next) => {
    req.is_authenticated = false;

    if (req.header("auth") || req.header("Authorization")) {
      const jwtToken = req.header("auth") || req.header("Authorization")?.split(" ")[1] || "";
      try {
        jwt.verify(jwtToken, JWT_SECRET, async function (err, decoded) {
          if (err || !decoded || !decoded.email) {
            console.log(err?.name);
            if (err?.name == "TokenExpiredError")
              return res.status(440).send("JWT authentication token expired. Please log in again");

            return res.status(401).send("Invalid authentication");
          }

          UserService.getUserById(decoded.id).then(
            (result) => {
              if (!result) {
                return res
                  .status(401)
                  .send("JWT authentication token invalid. Please log in again");
              }
              req.is_authenticated = true;
              req.user_id = decoded.id;
              req.email = result!.email;
              req.is_administrator = result!.isAdministrator;
              req.is_verified = result!.isEmailVerified;
              next();
            },
            (err) => {
              console.error(err);
              return res.status(500).send("Error during authentication: Failed to fetch user");
            }
          );
        });
      } catch (e) {
        return console.error(e); // abort on error;
      }
    } else {
      return res.status(403).send("JWT missing.");
    }
  });

  return router;
}
