import REST_Host from "./_rest_host";
import express, {Router} from "express";
import register_route from "./routes/register";
import verify_email from "./routes/verifyEmail";
import resend_verification_link from "./routes/resendVerificationLink";
import approve_user from "./routes/approveUser";

// setup the websocket-server on top of the http_server
export function express_routes(this:REST_Host) : Router {
    let router = express.Router();

    // TODO handle our REST routes
    this.expressApp.use(register_route());
    this.expressApp.use(verify_email());
    this.expressApp.use(resend_verification_link());
    this.expressApp.use(approve_user());
    this.expressApp.use((req, res) => res.status(404).send("Not found."));

    return router;
}