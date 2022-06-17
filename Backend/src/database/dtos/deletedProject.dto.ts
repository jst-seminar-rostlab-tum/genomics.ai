import { ObjectId } from "mongoose";

export interface AddDeletedProjectDTO {
  _id: ObjectId;
  owner: ObjectId;
  teamId: ObjectId;
  name: string;

  atlasId: ObjectId;
  modelId: ObjectId;

  uploadId?: string;
  location?: string;
  fileName?: string;
  fileSize?: number;
  uploadDate?: Date;

  status?: string;
  resultName?: string;
  resultSize?: number;

  deletedAt: Date;
}