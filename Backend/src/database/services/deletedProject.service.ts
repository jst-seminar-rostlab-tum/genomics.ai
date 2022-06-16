import { ObjectId } from "mongoose";
import { AddDeletedProjectDTO } from "../dtos/deletedProject.dto";
import {
  deletedProjectModel,
  DELETED_PROJECT_LIFETIME_DAYS,
  IDeletedProject,
} from "../models/deleted_projects";

export default class DeletedProjectService {
  /**
   *  Adds given project to the database.
   *
   *  @param    project
   *  @returns  projectAdded - the added project
   */
  static async addDeletedProject(
    project: AddDeletedProjectDTO
  ): Promise<IDeletedProject & { _id: ObjectId }> {
    return await deletedProjectModel.create(project);
  }

  /**
   * Search for the project with the given id
   * @param id project id
   * @returns the project
   */
  static async getDeletedProjectById(
    id: ObjectId | string
  ): Promise<(IDeletedProject & { _id: ObjectId }) | null> {
    return await deletedProjectModel.findById(id);
  }

  /**
   * Delete the project with the given id
   * @param id project id
   * @returns the deleted project
   */
  static async deleteDeletedProjectById(
    id: ObjectId | string
  ): Promise<(IDeletedProject & { _id: ObjectId }) | null> {
    return await deletedProjectModel.findByIdAndDelete(id);
  }

  /**
   *  List projects with the given owner and sort according to the two sort parameters
   *
   *  @param   owner - userId
   *  @param   sortUploadDate - sorting by upload date
   *  @param   sortDeletionDate - sorting by deletion date
   *  @returns project or null
   */
  static async getDeletedProjectsByOwner(
    user_id: ObjectId,
    sortUploadDate: number = 0,
    sortDeletionDate: number = 0
  ): Promise<(IDeletedProject & { _id: ObjectId })[]> {
    let sorting: any = {};
    if (sortUploadDate != 0) {
      sorting.uploadDate = Math.sign(sortUploadDate);
    }
    if (sortDeletionDate != 0) {
      sorting.deletionDate = Math.sign(sortDeletionDate);
    }
    //Return only projects which are less than LIFETIME_DAYS old => calculate the earliest deletion date after which projects should be returned
    let deletedAfter = new Date();
    deletedAfter.setTime(deletedAfter.getTime() - DELETED_PROJECT_LIFETIME_DAYS * 86_400_000);
    return await deletedProjectModel
      .find({ owner: user_id, deletedAt: { $gt: deletedAfter } })
      .sort(sorting);
  }

  static async getProjectsOverLifetime(): Promise<(IDeletedProject & { _id: ObjectId })[]> {
    let deletedBefore = new Date();
    deletedBefore.setTime(deletedBefore.getTime() - DELETED_PROJECT_LIFETIME_DAYS * 86_400_000);
    return await deletedProjectModel.find({ deletedAt: { $lt: deletedBefore } });
  }

  static async deleteProjectsByIds(projects: ObjectId[]) {
    return await deletedProjectModel.deleteMany({ _id: { $in: projects } });
  }
}
