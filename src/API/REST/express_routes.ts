import REST_Host from "./_rest_host";
import express, {Router} from "express";
import register_route from "./routes/register";
import auth_route from "./routes/auth";
import hello_route from "./routes/hello";
import hello_auth_route from "./routes/hello_auth";
import update_profile_route from "./routes/update_profile";
import get_profile_route from "./routes/get_profile";
import get_job_route from "./routes/get_job";
import get_jobs_route from "./routes/get_jobs";
import get_unauthorized_users_route from "./routes/get_unauthorized_users";
import authorize_user_route from "./routes/authorize_user";
import verify_email_route from "./routes/verify_email";
 
import resend_verification_link from "./routes/resend_verification_link";
// setup the websocket-server on top of the http_server
export function express_routes(this:REST_Host) : Router {
    let router = express.Router();

    // unauthenticated routes
    this.expressApp.use(auth_route());
    this.expressApp.use(register_route());
    this.expressApp.use(resend_verification_link());

    // authenticated routes
    this.expressApp.use(update_profile_route());
    this.expressApp.use(get_profile_route());
    this.expressApp.use(get_job_route());
    this.expressApp.use(get_jobs_route());

    // administrator routes
    this.expressApp.use(get_unauthorized_users_route());
    this.expressApp.use(authorize_user_route());

    // debugging / testing routes
    this.expressApp.use(hello_route());
    this.expressApp.use(hello_auth_route());

    this.expressApp.use(/^.*_ah.*$/, (req, res)=>res.status(200).send()) // always tell google everything is fine
    this.expressApp.use((req, res) => res.status(404).send("Not found."));

    return router;
}
