import { IDemo, demoModel } from "../models/demo";
import { ObjectId } from 'mongoose';

export default class demoService {
    /**
     * Get all the available demos.
     * 
     * @returns demo array
     */
    static async getAllDemos(): Promise<(IDemo & { _id: ObjectId})[]> {
        return demoModel.find().exec();
    }
}