import {ObjectId} from "mongoose";

/**
 *  Stores the raw data needed to add a project.
 */
export interface AddProjectDTO {
    title: string,
    description: string,
    visibility: string,
    adminIds: Array<ObjectId>,
}
