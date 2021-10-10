import express, {Router} from "express";

export default function register_route(): Router {
    let router = express.Router();

    router
        .post("/register")
        .use(((req, res) => {
            res.send("User has been registered");
        }))

    return router;
}