import { IInstitution, institutionModel } from "../models/institution";
import { AddInstitutionDTO, UpdateInstitutionDTO } from "../dtos/institution.dto";
import { ObjectId } from "mongoose";

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
        let institutionAdded: (IInstitution | undefined) = undefined;
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
        let updatedInstitution: (IInstitution | undefined) = undefined;

        const institution = await institutionModel.findOne({ _id: institutionId })

        if (institution) {

            institution.invitedMemberIds = [...institution.invitedMemberIds, userId];
            updatedInstitution = await institution.save();
            return updatedInstitution;
        } else {
            return undefined;
        }
    }

    /**
     *  Make a person admin of an institution.
     *
     *  @param    institutionId - the institution that is being modified 
     *  @param    userId - the user that is being added as admin
     *  @returns  institutionUpdated - institution with new admin if updated without error
     */
    static async makeUserAnAdminOfInstitution(institutionId: ObjectId, userId: ObjectId): Promise<IInstitution | undefined> {
        let updatedInstitution: (IInstitution | undefined) = undefined;

        const institution = await institutionModel.findOne({ _id: institutionId })

        if (institution) {

            institution.adminIds = [...institution.adminIds, userId];
            institution.memberIds = institution.memberIds.filter(id => String(id) !== String(userId))
            updatedInstitution = await institution.save();
            return updatedInstitution;
        } else {
            return undefined;
        }
    }

    /**
    *  Make a user member of an institution.
    *
    *  @param    institutionId - the institution that is being modified 
    *  @param    userId - the user that is being added as admin
    *  @returns  institutionUpdated - institution with new member if updated without error
    */
    static async makeUserMemberOfInstitution(institutionId: ObjectId, userId: ObjectId): Promise<IInstitution | undefined> {
        let updatedInstitution: (IInstitution | undefined) = undefined;

        const institution = await institutionModel.findOne({ _id: institutionId })

        if (institution) {

            institution.memberIds = [...institution.memberIds, userId];
            institution.invitedMemberIds = institution.invitedMemberIds.filter(id => String(id) !== String(userId))
            updatedInstitution = await institution.save();
            return updatedInstitution;
        } else {
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
    static async findMemeberOrInvitedById(user_id: (ObjectId | string), institution_id: (ObjectId | string)):
        Promise<(IInstitution & { _id: any; }) | undefined> {
        const result = await institutionModel.findOne({
            _id: institution_id,
            $or: [{
                invitedMemberIds: { $elemMatch: { $eq: user_id } },
            },
            {
                memberIds: { $elemMatch: { $eq: user_id } }
            }]

        })
        if (result) {
            return result;
        } else {
            return undefined;
        }
    }

    /**
     *  Search for member of institution by id if they exists.
     *
     *  @param   user_id - the user id to search for
     *  @param   institution_id - the institution id to search for
     *  @returns institution - if user is member of the institution
     */
    static async findMemeberById(user_id: (ObjectId | string), institution_id: (ObjectId | string)):
        Promise<(IInstitution & { _id: any; }) | undefined> {
        const result = await institutionModel.findOne({
            _id: institution_id,
            $or: [{
                memberIds: { $elemMatch: { $eq: user_id } }
            }]

        })

        if (result) {
            return result;
        } else {
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
        Promise<(IInstitution & { _id: ObjectId } | null)> {
        return await institutionModel.findOne({ name });
    }

    /**
     *  Search for an institution with the given id and return if found.
     *
     *  @param   id - the institution id to search for
     *  @returns institution - matching institution for id or null
     */
    static async getInstitutionById(institutionId: (ObjectId | string)):
        Promise<(IInstitution & { _id: ObjectId } | null)> {
        return await institutionModel.findById(institutionId).exec();
    }

    /**
     *  Updates the given institution corresponding to the id with the
     *  update_object.
     *
     *  @param institution_id
     *  @param update_object - includes fields to be updated
     */
    static async updateInstitution(institution_id: (ObjectId | string), update_object: UpdateInstitutionDTO) {
        await institutionModel.updateOne({ _id: institution_id }, update_object);
    }

    static async unsetProfilePicture(institution_id: ObjectId | string): Promise<string | null | undefined> {
        let old = await institutionModel.findByIdAndUpdate(institution_id, { $unset: { profilePictureURL: "" } });
        return old?.profilePictureURL;
    }
    static async unsetBackgroundPicture(institution_id: ObjectId | string): Promise<string | null | undefined> {
        let old = await institutionModel.findByIdAndUpdate(institution_id, { $unset: { backgroundPictureURL: "" } });
        return old?.backgroundPictureURL;
    }

    static async filterInstitutions(query: any):
        Promise<(IInstitution[] | null)> {

        var keyword: object,
            sortBy: any;

        query.hasOwnProperty('keyword') ? keyword = { name: query.keyword } : keyword = {};

        if (query.hasOwnProperty('sortBy')) {
            let sortProperty = query.sortBy;
            sortBy = { sortProperty: 1 }
        } else
            sortBy = {};

        return await institutionModel.find(keyword).sort(sortBy).exec();
    }


    /**
     *  Remove the given userId from the given institution.
     *
     *  @param   institutionId
     *  @param   userId
     *  @returns updateDocument
     */
     static async removeMemberFromTeam(institutionId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await institutionModel.updateOne(
            { _id: institutionId },
            { 
                $pull: { memberIds: userId, adminIds: userId}
            }
        );
    }
}
