import {IProject, projectModel} from "../models/project";
import {AddProjectDTO} from "../dtos/project.dto";
import {ObjectId} from "mongoose";

/**
 *  @class ProjectService
 *
 *  Provides useful methods to access the database and modify projects,
 *  which can be used by the route-controllers.
 */
export default class ProjectService {
    /**
     *  Adds given project to the database.
     *
     *  @param    project
     *  @returns  projectAdded - the added project
     */
    static async addProject(project: AddProjectDTO): Promise<IProject> {
        let projectAdded : (IProject | undefined) = undefined;
        projectAdded = await projectModel.create(project);
        return projectAdded;
    }

    /**
     *  Search for a project with the given title and return if found.
     *
     *  @param   title
     *  @returns project or null
     */
    static async getProjectByTitle(title: string):
      Promise<( IProject & { _id: ObjectId } | null )> {
        return await projectModel.findOne({title});
    }
}
