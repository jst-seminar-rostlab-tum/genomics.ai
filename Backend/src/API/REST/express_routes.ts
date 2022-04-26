import REST_Host from "./_rest_host";
import express, {Router} from "express";
import register_route from "./routes/register";
import auth_route from "./routes/auth";
import hello_route from "./routes/hello";
import hello_auth_route from "./routes/hello_auth";
import update_profile_route from "./routes/update_profile";
import initiate_processing_route from "./routes/initiate_processing";
import abort_processing_route from "./routes/abort_processing";

import get_profile_route from "./routes/get_profile";
import get_project_route from "./routes/get_project";
import get_projects_route from "./routes/get_projects";
import get_unauthorized_users_route from "./routes/get_unauthorized_users";
import authorize_user_route from "./routes/authorize_user";
import verify_email_route from "./routes/verify_email";
import password_reset_route from "./routes/password_reset";
 
import resend_verification_link from "./routes/resend_verification_link";
import upload_complete_upload_route from "./routes/file_upload/complete_upload";
import upload_start_upload_route from "./routes/file_upload/start_upload";
import upload_get_upload_url_route from "./routes/file_upload/get_upload_url";
import download_results_route from "./routes/file_download/results";

import {test_institution, create_institution, invite_to_institution} from "./routes/institution/institutionRouter";
import {create_team, invite_person_to_a_team, add_user_to_admin, join_member, 
    add_team_to_institution} from "./routes/team/teamRouter";

// setup the websocket-server on top of the http_server
export function express_routes(this: REST_Host): Router {
    let router = express.Router();

    // unauthenticated routes
    this.expressApp.use(auth_route());
    this.expressApp.use(register_route());
    this.expressApp.use(resend_verification_link());
    this.expressApp.use(verify_email_route());
    this.expressApp.use(password_reset_route());

    // authenticated routes
    this.expressApp.use(update_profile_route());
    this.expressApp.use(get_profile_route());
    this.expressApp.use(get_project_route());
    this.expressApp.use(get_projects_route());

    // administrator routes
    this.expressApp.use(get_unauthorized_users_route());
    this.expressApp.use(authorize_user_route());

    // debugging / testing routes
    this.expressApp.use(hello_route());
    this.expressApp.use(hello_auth_route());
    this.expressApp.use(initiate_processing_route());
    this.expressApp.use(abort_processing_route());

    // upload routes
    this.expressApp.use(upload_get_upload_url_route());
    this.expressApp.use(upload_start_upload_route());
    this.expressApp.use(upload_complete_upload_route());

    // download routes
    this.expressApp.use(download_results_route());

    // institution routes
    this.expressApp.use(create_institution());
    this.expressApp.use(invite_to_institution());
    this.expressApp.use(test_institution());

    // team routes
    this.expressApp.use(create_team());
    this.expressApp.use(invite_person_to_a_team());
    this.expressApp.use(add_user_to_admin());
    this.expressApp.use(join_member());
    this.expressApp.use(add_team_to_institution());

    this.expressApp.use(/^.*_ah.*$/, (req, res) => res.status(200).send()) // always tell google everything is fine
    this.expressApp.use((req, res) => res.status(404).send("Not found."));

    return router;
}
