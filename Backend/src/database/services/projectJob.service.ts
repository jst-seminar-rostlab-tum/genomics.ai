import {IProjectJob, projectJobModel} from "../models/projectJob";
import {AddProjectJobDTO, UpdateProjectJobDTO} from "../dtos/projectJob.dto";
import {ObjectId} from "mongoose";

/**
 *  @class ProjectJobService
 *
 *  Provides useful methods to access the database and modify project jobs,
 *  which can be used by the route-controllers.
 */
export default class ProjectJobService {
    /**
     *  Search for a project job with the given project id and return if found.
     *
     *  @param   project_id - the project id to search for
     *  @returns project job - matched project job to project_id or null
     */
    static async getProjectJobById(project_id: (ObjectId | string)):
      Promise<( IProjectJob & { _id: ObjectId } | null )> {
        return await projectJobModel.findById(project_id).exec();
    }

    /**
     *  Search for a project job with the given uploadId and project owner (userId - optional)
     *  and return if found.
     *
     *  @param   uploadId
     *  @param   owner? - userId
     *  @returns project job or null
     */
    static async getProjectJobByUploadId(uploadId: string, owner?: ObjectId):
      Promise<( IProjectJob & { _id: ObjectId } | null )> {
        return typeof(owner) === 'undefined' ?
          await projectJobModel.findOne({ uploadId }).exec() :
          await projectJobModel.findOne({ uploadId, owner }).exec();
    }

    /**
     *  Search for a project job with the given project owner and sort
     *  in order of the given sort parameter.
     *
     *  @param   owner - userId
     *  @param   sort - order of the sort
     *  @returns project jobs or null
     */
    static async getProjectJobByOwner(user_id: ObjectId, sort: number = 0):
      Promise<( IProjectJob & { _id: ObjectId } )[]> {
        return await projectJobModel.find({owner: user_id}).sort({uploadDate: sort});
    }

    /**
     *  Adds given project job to the database.
     *
     *  @param    project job
     *  @returns  projectJobAdded - the added project job
     */
    static async addProjectJob(projectJob: AddProjectJobDTO): Promise<IProjectJob> {
        let projectJobAdded : (IProjectJob | undefined) = undefined;
        projectJobAdded = await projectJobModel.create(projectJob);
        return projectJobAdded;
    }

    /**
     *  Updates the given project job corresponding to the uploadId with the
     *  update_object.
     *
     *  @param uploadId
     *  @param update_object - includes fields to be updated
     */
    static async updateProjectJob(uploadId: string, update_object: UpdateProjectJobDTO) {
        await projectJobModel.updateOne({uploadId}, update_object).exec();
    }

    /**
     *  Updates the upload id of the given project job.
     *
     *  @param _id
     *  @param uploadId
     */
    static async updateUploadId(_id: ObjectId, uploadId: string) {
        await projectJobModel.updateOne({_id}, {uploadId}).exec();
    }
}
