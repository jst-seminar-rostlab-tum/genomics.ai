import { ObjectId } from "mongoose";

/**
 *  Stores the raw data to update a project.
 */
export interface UpdateProjectDTO {
  fileSize?: number;
  status?: string;
  location?: string;
}

/**
 * Stores the raw data needed to create a project.
 */
export interface AddProjectDTO {
    owner: ObjectId,
    name: string,
    fileName: string,
    uploadDate: Date,
    status: string,
    modelId: string,
    atlasId: string,
}
