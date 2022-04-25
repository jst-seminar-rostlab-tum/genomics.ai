import {ObjectId} from "mongoose";

/**
 *  Stores the raw data needed to add a team.
 */
export interface AddTeamDTO {
    title: string,
    description: string,
    visibility: string,
    adminIds: Array<ObjectId>,
}
