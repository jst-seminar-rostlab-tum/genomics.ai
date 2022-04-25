import {IModel, modelModel} from "../models/model";
import {ObjectId} from "mongoose";

export default class ModelService {
    /**
     *  Search for a team with the given team id and return if found.
     *
     *  @param   modelId
     *  @returns model - matched model to modelId or null
     */
    static async getModelById(modelId: (ObjectId | string)):
      Promise<( IModel & { _id: ObjectId } | null )> {
        return await modelModel.findById(modelId).exec();
    }
}
