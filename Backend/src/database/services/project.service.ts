import { IProject, projectModel } from "../models/project";
import { AddProjectDTO, UpdateProjectDTO } from "../dtos/project.dto";
import { ObjectId } from "mongoose";

/**
 *  @class ProjectService
 *
 *  Provides useful methods to access the database and modify projects,
 *  which can be used by the route-controllers.
 */
export default class ProjectService {
  /**
   *  Search for a project with the given team id and return if found.
   *
   *  @param   project_id - the team id to search for
   *  @returns project - matched project to project_id or null
   */
  static async getProjectById(
    project_id: ObjectId | string
  ): Promise<(IProject & { _id: ObjectId }) | null> {
    return await projectModel.findById(project_id).exec();
  }

  /**
   * Deletes a project with the given project id and return the deleted document.
   * 
   * @param project_id - the project id
   * @returns project - project or null if not found
   */
  static async deleteProjectById(
    project_id: ObjectId | string
  ): Promise<(IProject & { _id: ObjectId}) | null> {
    return await projectModel.findByIdAndRemove(project_id).exec();
  }

  /**
   *  Search for a project with the given uploadId and team owner (userId - optional)
   *  and return if found.
   *
   *  @param   uploadId
   *  @param   owner? - userId
   *  @returns project or null
   */
  static async getProjectByUploadId(
    uploadId: string,
    owner?: ObjectId
  ): Promise<(IProject & { _id: ObjectId }) | null> {
    return typeof owner === "undefined"
      ? await projectModel.findOne({ uploadId }).exec()
      : await projectModel.findOne({ uploadId, owner }).exec();
  }

  /**
   *  Search for a project with the given team owner and sort
   *  in order of the given sort parameter.
   *
   *  @param   owner - userId
   *  @param   sort - order of the sort
   *  @returns project or null
   */
  static async getProjectByOwner(
    user_id: ObjectId,
    sort: number = 0
  ): Promise<(IProject & { _id: ObjectId })[]> {
    return await projectModel.find({ owner: user_id }).sort({ uploadDate: sort });
  }

  /**
   *  Returns projects of the given team ids,
   *  where the teams of the projects are in
   *  the given team list.
   *
   *  @param   teamIds as ObjectId array
   *  @returns projects or null
   */
  static async getProjectsOfTeams(
    teamIds: ObjectId[]
  ): Promise<(IProject & { _id: ObjectId })[]> {
    return await projectModel.find({ teamId: { $in: teamIds } }).exec();
  }

  /**
   *  Returns projects of the given user ids,
   *  where the owners of the projects are in
   *  the given user list.
   *
   *  @param   userIds as ObjectId array
   *  @returns projects or null
   */
  static async getProjectsOfUsers(
    userIds: ObjectId[]
  ): Promise<(IProject & { _id: ObjectId })[]> {
    return await projectModel.find({ owner: { $in: userIds } }).exec();
  }

  /**
   *  Adds given project to the database.
   *
   *  @param    project
   *  @returns  projectAdded - the added project
   */
  static async addProject(project: AddProjectDTO): Promise<IProject> {
    let projectAdded: IProject | undefined = undefined;
    projectAdded = await projectModel.create(project);
    return projectAdded;
  }

  /**
   *  Updates the given project corresponding to the uploadId with the
   *  update_object.
   *
   *  @param uploadId
   *  @param update_object - includes fields to be updated
   */
  static async updateProjectByUploadId(uploadId: string, update_object: UpdateProjectDTO) {
    await projectModel.updateOne({ uploadId }, update_object).exec();
  }

  /**
   *  Updates the project with id with update_object.
   *
   *  @param uploadId
   *  @param update_object - includes fields to be updated
   */
  static async updateProjectById(id: ObjectId | string, update_object: UpdateProjectDTO) {
    await projectModel.findByIdAndUpdate(id, update_object);
  }

  /**
   *  Updates the upload id of the given project.
   *
   *  @param _id
   *  @param uploadId
   */
  static async updateUploadId(_id: ObjectId, uploadId: string) {
    await projectModel.updateOne({ _id }, { uploadId }).exec();
  }

  /**
   *  Add the given userId to the admin list and removes he/she from the memberIds, of the given project.
   *
   *  @param   projectId
   *  @param   userId
   *  @returns updateDocument
   */
  static async addAdminToProject(
    projectId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await projectModel.updateOne(
      { _id: projectId },
      {
        $addToSet: { adminIds: userId },
        $pull: { memberIds: userId },
      }
    );
  }

  /**
   *  Add the given userId into the given project.
   *
   *  @param   projectId
   *  @param   userId
   *  @returns updateDocument
   */
  static async addNewMemberIntoProject(
    projectId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await projectModel.updateOne({ _id: projectId }, { $addToSet: { memberIds: userId } });
  }

  /**
   *  Set the team of a project to the given teamId.
   *
   *  @param   projectId
   *  @param   teamId
   *  @returns updateDocument
   */
  static async setTeamOfProject(
    projectId: ObjectId | string,
    teamId: ObjectId | string
  ): Promise<any> {
    return await projectModel.updateOne({ _id: projectId }, { $set: { teamId: teamId } });
  }

  static async getProjects(queryParams: any): Promise<IProject[] | null> {
    var keyword: object, sortBy: any;

    queryParams.hasOwnProperty("keyword")
      ? (keyword = { name : { $regex : "^" +  queryParams.keyword, $options : 'i'}})
      : (keyword = {});

    if (queryParams.hasOwnProperty("sortBy")) {
      let sortProperty = queryParams.sortBy;
      sortBy = { sortProperty: 1 };
    } else sortBy = {};

    return await projectModel.find(keyword).sort(sortBy);
  }
}
