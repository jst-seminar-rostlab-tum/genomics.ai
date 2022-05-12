import { ObjectId } from "mongoose";

/**
 *  Stores the raw data needed to create a project update token.
 */
export interface AddProjectUpdateTokenDTO {
  _projectId: ObjectId;
}