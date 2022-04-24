import {ObjectId} from "mongoose";

/**
 *  Stores the raw data needed to add an institution.
 */
export interface AddInstitutionDTO {
    name:string,
    country:string,
    profilePictureURL: string,
    backgroundPictureURL: string,
    adminIds: Array<ObjectId>,
}

export interface UpdateInstitutionDTO {
    name?:string,
    country?: string,
    profilePictureURL?:string,
    backgroundPictureUrl?: string,
}
