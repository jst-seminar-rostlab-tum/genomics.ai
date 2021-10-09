import express from "express";

export default function auth_route(){
    let router = express.Router();

    router.post("/auth", (req, res, next) => {

        res.send("Hello World")
    })

    return router;
}