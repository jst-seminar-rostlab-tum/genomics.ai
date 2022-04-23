import {IInstitution, institutionModel} from "../models/institution";
import {AddInstitutionDTO} from "../dtos/institution.dto";
import {ObjectId} from "mongoose";
import { IUser } from "../models/user";

/**
 *  @class InstitutionService
 *
 *  Provides useful methods to access the database and modify institutions,
 *  which can be used by the route-controllers.
 */
export default class InstitutionService {
    /**
     *  Adds given institution to the database.
     *
     *  @param    institution
     *  @returns  institutionAdded - the added institution
     */
    static async addInstitution(institution: AddInstitutionDTO): Promise<IInstitution> {
        let institutionAdded : (IInstitution | undefined) = undefined;
        institutionAdded = await institutionModel.create(institution);
        return institutionAdded;
    }

    /**
     *  Invite a person to an institution.
     *
     *  @param    institutionId
     *  @param    userId
     *  @returns  institutionUpdated 
     */
     static async inviteToInstitution(institutionId: ObjectId, userId: ObjectId): Promise<IInstitution | undefined> {
        let updatedInstitution : (IInstitution | undefined) = undefined;

        const institution = await institutionModel.findOne({_id: institutionId})
                
        if (institution) {

            institution.invitedMemberIds = [...institution.invitedMemberIds, userId];
            updatedInstitution = await institution.save();
            return updatedInstitution;
        }else {
            return undefined;
        }

    }

    /**
     *  Search for invited member or member of institution by id if they exists.
     *
     *  @param   user_id - the user id to search for
     *  @param   institution_id - the institution id to search for
     *  @returns institution - if user is member of the institution
     */
     static async findMemeberById(user_id: (ObjectId | string), institution_id: (ObjectId | string)):
     Promise<(IInstitution & { _id: any; }) | undefined>  {
        const result = await institutionModel.findOne({
            _id: institution_id, 
            $or: [{
                invitedMemberIds: { $elemMatch: {$eq: user_id} },
            },
            {
                memberIds: { $elemMatch: {$eq: user_id} }
            }]
            
        })
        console.log(result)
        if (result) {
            return result;
        }else {
            return undefined;
        }
    }

    /**
     *  Search for an institution with the given name and return if found.
     *
     *  @param   name
     *  @returns institution or null
     */
    static async getInstitutionByName(name: string):
      Promise<( IInstitution & { _id: ObjectId } | null )> {
        return await institutionModel.findOne({name});
    }
}
