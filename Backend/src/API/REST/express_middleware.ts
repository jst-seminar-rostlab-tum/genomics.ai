import REST_Host from "./_rest_host";
import express, { Router } from "express";
import rateLimit from "express-rate-limit";
import bodyParser from "body-parser";
import cors from "cors";

// express.js middleware before requests hit route-handlers
export function express_middleware(this: REST_Host): Router {
  let router = express.Router();

  router.use(cors());
  router.use(
    rateLimit({
      windowMs: 60 * 1000, // 1min
      max: 1000, // limit each IP to 1000 requests per windowMs
      message: "Limit exceeded. Try again in an hour",
    })
  );

  router.use(bodyParser.json());
  router.use(bodyParser.urlencoded({ extended: true }));
  return router;
}
export function express_routes_middleware(this: REST_Host): Router {
  let router = express.Router();
  router.use(
    ["/user-avatar", "/institutions/:id/profilepicture"],
    bodyParser.raw({
      limit: "5MB",
      type: ["image/png", "image/jpeg"],
    })
  );
  router.use(
    "/institutions/:id/backgroundpicture",
    bodyParser.raw({
      limit: "8MB",
      type: ["image/png", "image/jpeg"],
    })
  );

  return router;
}
