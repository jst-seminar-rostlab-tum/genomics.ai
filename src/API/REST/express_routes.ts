import REST_Host from "./_rest_host";
import express, {Router} from "express";
import register_route from "./routes/register";
import auth_route from "./routes/auth";
import hello_route from "./routes/hello";
import hello_auth_route from "./routes/hello_auth";
import update_profile_route from "./routes/update_profile";
import initiate_processing_route from "./routes/initiate_processing";
import abort_processing_route from "./routes/abort_processing";

// setup the websocket-server on top of the http_server
export function express_routes(this:REST_Host) : Router {
    let router = express.Router();

    // unauthenticated routes
    this.expressApp.use(auth_route());
    this.expressApp.use(register_route());

    // authenticated routes
    this.expressApp.use(update_profile_route());

    // debugging / testing routes
    this.expressApp.use(hello_route());
    this.expressApp.use(hello_auth_route());
    this.expressApp.use(initiate_processing_route());
    this.expressApp.use(abort_processing_route());

    this.expressApp.use(/^.*_ah.*$/, (req, res)=>res.status(200).send()) // always tell google everything is fine
    this.expressApp.use((req, res) => res.status(404).send("Not found."));

    return router;
}