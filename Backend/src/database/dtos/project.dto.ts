import {ObjectId} from "mongoose";

/**
 *  Stores the raw data to update a team job.
 */
export interface UpdateProjectDTO {
    fileSize?: number;
    status?: string;
    location?: string;
}

/**
 * Stores the raw data needed to create a team job.
 */
export interface AddProjectDTO {
    owner: ObjectId,
    fileName: string,
    uploadDate: Date,
    status: string
}
