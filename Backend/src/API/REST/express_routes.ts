import REST_Host from "./_rest_host";
import express, { Router } from "express";
import register_route from "./routes/register";
import auth_route from "./routes/auth";
import hello_route from "./routes/hello";
import hello_auth_route from "./routes/hello_auth";
import update_profile_route from "./routes/update_profile";
import initiate_processing_route from "./routes/initiate_processing";
import abort_processing_route from "./routes/abort_processing";

import get_profile_route from "./routes/get_profile";
//import get_project_route from "./routes/get_project";
//import get_projects_route from "./routes/get_projects";
import get_unauthorized_users_route from "./routes/get_unauthorized_users";
import authorize_user_route from "./routes/authorize_user";
import verify_email_route from "./routes/verify_email";
import password_reset_route from "./routes/password_reset";

import resend_verification_link from "./routes/resend_verification_link";
import upload_complete_upload_route from "./routes/file_upload/complete_upload";
import upload_start_upload_route from "./routes/file_upload/start_upload";
import upload_get_upload_url_route from "./routes/file_upload/get_upload_url";
import download_results_route from "./routes/file_download/results";
import upload_user_avatar_route from "./routes/upload_user_avatar";

import { get_teams_of_user, get_users } from "./routes/user/userRouter";
import { get_model, get_allModels } from "./routes/model/modelRouter";
import { get_atlas, get_atlas_visualization, get_allAtlases } from "./routes/atlas/atlasRouter";
import * as swaggerDocument from "../../swagger.json";
import * as swaggerUi from "swagger-ui-express";

import {
  create_institution,
  invite_to_institution,
  make_user_admin_of_institution,
  join_as_member_of_institution,
  get_institutions,
  get_institution,
  get_members_of_institution,
  get_teams_of_institution,
  get_projects_of_institution,
  get_users_institutions,
  disjoin_member_of_institution,
} from "./routes/institution/institutionRouter";

import {
  create_team,
  invite_person_to_a_team,
  add_user_to_admin,
  join_member,
  add_team_to_institution,
  remove_team_from_institution,
  add_project_to_team,
  get_teams,
  get_users_teams,
  disjoin_member,
  get_team,
  update_team,
} from "./routes/team/teamRouter";

import {
  get_projects,
  get_userProjects,
  get_project_by_id,
  get_users_projects,
} from "./routes/project/projectRouter";

import {
  upload_institution_backgroundpicture_route,
  upload_institution_profilepicture_route,
} from "./routes/upload_institution_pictures";

import {
  reset_institution_backgroundpicture_route,
  reset_institution_profilepicture_route,
} from "./routes/reset_institution_pictures";
import reset_user_avatar_route from "./routes/reset_user_avatar";

// setup the websocket-server on top of the http_server
export function express_routes(this: REST_Host): Router {
  let router = express.Router();

  this.expressApp.use("/api-docs", swaggerUi.serve, swaggerUi.setup(swaggerDocument));

  // unauthenticated routes
  this.expressApp.use(auth_route());
  this.expressApp.use(register_route());
  this.expressApp.use(resend_verification_link());
  this.expressApp.use(verify_email_route());
  this.expressApp.use(password_reset_route());

  // authenticated routes
  this.expressApp.use(update_profile_route());
  this.expressApp.use(get_profile_route());
  //this.expressApp.use(get_project_route());  -- depreciated, we use get_project_by_id
  //this.expressApp.use(get_projects_route()); -- depreciated, we use get_userProjects
  this.expressApp.use(upload_user_avatar_route());
  this.expressApp.use(reset_user_avatar_route());

  // administrator routes
  this.expressApp.use(get_unauthorized_users_route());
  this.expressApp.use(authorize_user_route());

  // debugging / testing routes
  this.expressApp.use(hello_route());
  this.expressApp.use(hello_auth_route());
  //this.expressApp.use(initiate_processing_route());
  //this.expressApp.use(abort_processing_route());

  // team routes
  this.expressApp.use(create_team());
  this.expressApp.use(invite_person_to_a_team());
  this.expressApp.use(add_project_to_team());
  this.expressApp.use(add_user_to_admin());
  this.expressApp.use(join_member());
  this.expressApp.use(add_team_to_institution());
  this.expressApp.use(remove_team_from_institution());
  this.expressApp.use(get_teams());
  this.expressApp.use(get_users_teams());
  this.expressApp.use(disjoin_member());
  this.expressApp.use(get_team());
  this.expressApp.use(update_team());

  // user routes
  this.expressApp.use(get_teams_of_user());
  this.expressApp.use(get_users());

  // project routes
  this.expressApp.use(get_projects());
  this.expressApp.use(get_userProjects());
  this.expressApp.use(get_project_by_id());
  this.expressApp.use(get_users_projects());

  // model routes
  this.expressApp.use(get_model());
  this.expressApp.use(get_allModels());

  // atlas routes
  this.expressApp.use(get_atlas());
  this.expressApp.use(get_atlas_visualization());
  this.expressApp.use(get_allAtlases());

  // upload routes
  this.expressApp.use(upload_get_upload_url_route());
  this.expressApp.use(upload_start_upload_route());
  this.expressApp.use(upload_complete_upload_route());

  // download routes
  this.expressApp.use(download_results_route());

  // institution routes
  this.expressApp.use(create_institution());
  this.expressApp.use(invite_to_institution());
  this.expressApp.use(get_institution());
  this.expressApp.use(get_institutions());
  this.expressApp.use(upload_institution_profilepicture_route());
  this.expressApp.use(upload_institution_backgroundpicture_route());
  this.expressApp.use(reset_institution_profilepicture_route());
  this.expressApp.use(reset_institution_backgroundpicture_route());
  this.expressApp.use(make_user_admin_of_institution());
  this.expressApp.use(join_as_member_of_institution());
  this.expressApp.use(get_members_of_institution());
  this.expressApp.use(get_teams_of_institution());
  this.expressApp.use(get_projects_of_institution());
  this.expressApp.use(get_users_institutions());
  this.expressApp.use(disjoin_member_of_institution());

  this.expressApp.use(/^.*_ah.*$/, (req, res) => res.status(200).send()); // always tell google everything is fine
  this.expressApp.use((req, res) => res.status(404).send("Not found."));

  return router;
}
