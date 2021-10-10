import REST_Host from "./_rest_host";
import express, {Router} from "express";
import register_route from "./routes/register";

// setup the websocket-server on top of the http_server
export function express_routes(this:REST_Host) : Router {
    let router = express.Router();

    // TODO handle our REST routes
    this.expressApp.use(register_route());
    this.expressApp.use((req, res) => res.status(404).send("Not found."));

    return router;
}