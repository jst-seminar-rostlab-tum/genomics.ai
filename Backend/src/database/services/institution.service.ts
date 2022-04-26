import {IInstitution, institutionModel} from "../models/institution";
import {AddInstitutionDTO} from "../dtos/institution.dto";
import {ObjectId} from "mongoose";

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

    /**
     *  Search for a institution with the given institution id and return if found.
     *
     *  @param   institutionId
     *  @returns project - matched proejct to projectId or null
     */
    static async getInstitutionById(institutionId: (ObjectId | string)):
     Promise<( IInstitution & { _id: ObjectId } | null )> {
       return await institutionModel.findById(institutionId).exec();
    }

    /**
     *  Return the visibility status {"PRIVATE", "PUBLIC"}
     *  of the given institution if it exists, otherwise null.
     *
     *  @param  institutionId
     *  @return visibilityStatus - {"PRIVATE", "PUBLIC"}
     */
    static async getVisibility(institutionId: (ObjectId | string)):
      Promise<(string | null)> {
        const institution = await this.getInstitutionById(institutionId);
        if (!institution)
            return null;
        return institution.visibility;
    }

    /**
     *  Add the given teamId into the institution.
     *
     *  @param   teamId
     *  @param   institutionId
     *  @returns updateDocument
     */
     static async addNewMemberIntoTeam(teamId: (ObjectId | string), institutionId: (ObjectId | string)): Promise<any> {
        return await institutionModel.updateOne(
            { _id: institutionId },
            { $addToSet: { memberIds: teamId} }
        );
    }

    /**
     *  Add the given userId to the admin list and removes he/she from the memberIds, of the given team.
     *
     *  @param   teamId
     *  @param   institutionId
     *  @returns updateDocument
     */
     static async removeTeamFromInstitution(teamId: (ObjectId | string), institutionId: (ObjectId | string)): Promise<any> {
        return await institutionModel.updateOne(
            { _id: institutionId },
            {
                $pull: { memberIds: teamId }
            }
        );
    }

    /**
     *  Returns true if the given user is an admin of the given institution.
     *  The given institution should exist, otherwise the method returns false.
     *
     *  @param  userId
     *  @param  institutionId
     *  @return isAdmin
     */
    static async isAdmin(userId: (ObjectId | string), institutionId: (ObjectId | string)): Promise<boolean> {
        const institution = await this.getInstitutionById(institutionId);
        if (!institution)
          return false; /* institution does not exist */

        let isAdmin = false;
        var listAdmins = institution.adminIds.map(String);
        var userIdStr= String(userId);
        if (listAdmins.includes(userIdStr))
            isAdmin = true;
        return isAdmin;
    }

    /**
     *  Returns true if the given user is a member of the given institution.
     *  The given institution should exist, otherwise the method returns false.
     *
     *  @param  userId
     *  @param  institutionId
     *  @return isMember
     */
    static async isMember(userId: (ObjectId | string), institutionId: (ObjectId | string)): Promise<boolean> {
        const institution = await this.getInstitutionById(institutionId);
        if (!institution)
          return false; /* institution does not exist */

        let isMember = false;
        var listMembers = institution.memberIds.map(String);

        var userIdStr = String(userId);
        if (listMembers.includes(userIdStr))
            isMember = true;
        return isMember;
    }
}
