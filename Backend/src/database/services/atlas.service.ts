import {IAtlas, atlasModel} from "../models/atlas";
import {ObjectId} from "mongoose";

export default class AtlasService {
    /**
     *  Search for an atlas with the given atlas id and return if found.
     *
     *  @param   atlasId
     *  @returns atlas - matched atlas to atlasId or null
     */
    static async getAtlasById(atlasId: (ObjectId | string)):
      Promise<( IAtlas & { _id: ObjectId } | null )> {
        return await atlasModel.findById(atlasId).exec();
    }
}
