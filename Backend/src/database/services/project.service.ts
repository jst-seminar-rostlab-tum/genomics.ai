import {IProject, projectModel} from "../models/project";
import {AddProjectDTO, UpdateProjectDTO} from "../dtos/project.dto";
import {ObjectId} from "mongoose";

/**
 *  @class ProjectService
 *
 *  Provides useful methods to access the database and modify team jobs,
 *  which can be used by the route-controllers.
 */
export default class ProjectService {
    /**
     *  Search for a team job with the given team id and return if found.
     *
     *  @param   project_id - the team id to search for
     *  @returns project job - matched team job to project_id or null
     */
    static async getProjectById(project_id: (ObjectId | string)):
      Promise<( IProject & { _id: ObjectId } | null )> {
        return await projectModel.findById(project_id).exec();
    }

    /**
     *  Search for a team job with the given uploadId and team owner (userId - optional)
     *  and return if found.
     *
     *  @param   uploadId
     *  @param   owner? - userId
     *  @returns project job or null
     */
    static async getProjectByUploadId(uploadId: string, owner?: ObjectId):
      Promise<( IProject & { _id: ObjectId } | null )> {
        return typeof(owner) === 'undefined' ?
          await projectModel.findOne({ uploadId }).exec() :
          await projectModel.findOne({ uploadId, owner }).exec();
    }

    /**
     *  Search for a team job with the given team owner and sort
     *  in order of the given sort parameter.
     *
     *  @param   owner - userId
     *  @param   sort - order of the sort
     *  @returns project jobs or null
     */
    static async getProjectByOwner(user_id: ObjectId, sort: number = 0):
      Promise<( IProject & { _id: ObjectId } )[]> {
        return await projectModel.find({owner: user_id}).sort({uploadDate: sort});
    }

    /**
     *  Adds given team job to the database.
     *
     *  @param    team job
     *  @returns  projectJobAdded - the added team job
     */
    static async addProject(projectJob: AddProjectDTO): Promise<IProject> {
        let projectJobAdded : (IProject | undefined) = undefined;
        projectJobAdded = await projectModel.create(projectJob);
        return projectJobAdded;
    }

    /**
     *  Updates the given team job corresponding to the uploadId with the
     *  update_object.
     *
     *  @param uploadId
     *  @param update_object - includes fields to be updated
     */
    static async updateProject(uploadId: string, update_object: UpdateProjectDTO) {
        await projectModel.updateOne({uploadId}, update_object).exec();
    }

    /**
     *  Updates the upload id of the given team job.
     *
     *  @param _id
     *  @param uploadId
     */
    static async updateUploadId(_id: ObjectId, uploadId: string) {
        await projectModel.updateOne({_id}, {uploadId}).exec();
    }

    /**
     *  Returns the id of the owner of the given project if it exists,
     *  otherwise null.
     *
     *  @param   project_id
     *  @returns owner or null
     */
    static async getOwner(project_id: (ObjectId | string)):
      Promise<(ObjectId | null)> {
        const project = await ProjectService.getProjectById(project_id);
        if (!project)
            return null;
        return project.owner;
    }
}
