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

import { get_teams_of_user, get_users, get_user_by_id } from "./routes/user/userRouter";
import { get_model, get_allModels } from "./routes/model/modelRouter";
import { get_atlas, get_atlas_visualization, get_allAtlases } from "./routes/atlas/atlasRouter";

import * as swaggerUi from "swagger-ui-express";

import { loadSwaggerDocument } from "../../swagger/load-swagger";

import {
  create_institution,
  invite_to_institution,
  make_user_admin_of_institution,
  join_as_member_of_institution,
  get_institutions,
  get_institution,
  get_members_of_institution,
  remove_member_from_institution,
  remove_admin_role_for_institution_member,
  get_teams_of_institution,
  get_projects_of_institution,
  get_users_institutions,
  disjoin_member_of_institution,
  update_institution,
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
  get_members_of_team,
  remove_member_from_team,
  remove_admin_role_for_team_member,
  get_projects_of_team,
  remove_project_from_team,
} from "./routes/team/teamRouter";

import {
  get_projects,
  get_userProjects,
  get_project_by_id,
  get_users_projects,
  update_project_results,
  delete_project,
  get_deleted_projects,
  restore_deleted_project,
  cleanup_old_projects,
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
import { contact_us } from "./routes/contact/contactRoute";

import { get_allDemos } from "./routes/demo/demoRouter";

// setup the websocket-server on top of the http_server
export function express_routes(): Router {
  let router = express.Router();

  router.use("/api-docs", swaggerUi.serve, swaggerUi.setup(loadSwaggerDocument()));

  // unauthenticated routes
  router.use(auth_route());
  router.use(register_route());
  router.use(resend_verification_link());
  router.use(verify_email_route());
  router.use(password_reset_route());

  // authenticated routes
  router.use(update_profile_route());
  router.use(get_profile_route());
  //router.use(get_project_route());  -- depreciated, we use get_project_by_id
  //router.use(get_projects_route()); -- depreciated, we use get_userProjects
  router.use(upload_user_avatar_route());
  router.use(reset_user_avatar_route());

  // administrator routes
  router.use(get_unauthorized_users_route());
  router.use(authorize_user_route());

  // debugging / testing routes
  router.use(hello_route());
  router.use(hello_auth_route());
  //router.use(initiate_processing_route());
  //router.use(abort_processing_route());

  // team routes
  router.use(create_team());
  router.use(invite_person_to_a_team());
  router.use(add_project_to_team());
  router.use(add_user_to_admin());
  router.use(join_member());
  router.use(add_team_to_institution());
  router.use(remove_team_from_institution());
  router.use(get_teams());
  router.use(get_users_teams());
  router.use(disjoin_member());
  router.use(get_team());
  router.use(update_team());
  router.use(get_members_of_team());
  router.use(remove_member_from_team());
  router.use(remove_admin_role_for_team_member());
  router.use(get_projects_of_team());
  router.use(remove_project_from_team());

  // user routes
  router.use(get_teams_of_user());
  router.use(get_users());
  router.use(get_user_by_id());

  // project routes
  router.use(get_projects());
  router.use(get_userProjects());
  router.use(get_project_by_id());
  router.use(get_users_projects());
  router.use(delete_project());
  router.use(get_deleted_projects());
  router.use(restore_deleted_project());
  router.use(cleanup_old_projects());
  router.use(update_project_results());

  // model routes
  router.use(get_model());
  router.use(get_allModels());

  // atlas routes
  router.use(get_atlas());
  router.use(get_atlas_visualization());
  router.use(get_allAtlases());

  // demo routres
  router.use(get_allDemos());

  // upload routes
  router.use(upload_get_upload_url_route());
  router.use(upload_start_upload_route());
  router.use(upload_complete_upload_route());

  // download routes
  router.use(download_results_route());

  //contact routes
  router.use(contact_us());

  // institution routes
  router.use(create_institution());
  router.use(update_institution());
  router.use(invite_to_institution());
  router.use(get_institution());
  router.use(get_institutions());
  router.use(upload_institution_profilepicture_route());
  router.use(upload_institution_backgroundpicture_route());
  router.use(reset_institution_profilepicture_route());
  router.use(reset_institution_backgroundpicture_route());
  router.use(make_user_admin_of_institution());
  router.use(join_as_member_of_institution());
  router.use(get_members_of_institution());
  router.use(get_teams_of_institution());
  router.use(get_projects_of_institution());
  router.use(get_users_institutions());
  router.use(remove_member_from_institution());
  router.use(remove_admin_role_for_institution_member());
  router.use(disjoin_member_of_institution());

  router.use(/^.*_ah.*$/, (req, res) => res.status(200).send()); // always tell google everything is fine
  router.use((req, res) => res.status(404).send("Not found."));

  return router;
}
