import { ObjectId } from "mongoose";
import { AddProjectUpdateTokenDTO } from "../dtos/project_update_token.dto";
import { IProjectUpdateToken, projectUpdateTokenModel } from "../models/project_update_token";

/**
 *  @class ProjectUpdateTokenService
 *
 *  Provides useful methods to access the database and modify
 *  project-update-token, which can be used by the route-controllers.
 */
export default class ProjectUpdateTokenService {
  /**
   *  Creates project-update-token for the given userId.
   *
   *  @param    adddto - DTO with projectId field
   *  @returns  newly created project-update-token document
   */
  static async addToken(adddto: AddProjectUpdateTokenDTO): Promise<IProjectUpdateToken> {
    return await projectUpdateTokenModel.create({ _projectId: adddto._projectId });
  }

  /**
   *  Search project-update-token by the given token.
   *
   *  @param    token
   *  @returns  project-update-token or null
   */
  static async getTokenByToken(
    token: string
  ): Promise<(IProjectUpdateToken & { _id: ObjectId }) | null> {
    return await projectUpdateTokenModel.findOne({ token });
  }
}