import {ObjectId} from "mongoose";

/**
 *  Stores the raw data to update a project job.
 */
export interface UpdateProjectJobDTO {
    fileSize?: number;
    status?: string;
    location?: string;
}

/**
 * Stores the raw data needed to create a project job.
 */
export interface AddProjectJobDTO {
    owner: ObjectId,
    fileName: string,
    uploadDate: Date,
    status: string
}
