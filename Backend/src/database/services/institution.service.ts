import {IInstitution, institutionModel} from "../models/institution";
import {AddInstitutionDTO, UpdateInstitutionDTO} from "../dtos/institution.dto";
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
     *  Search for an institution with the given id and return if found.
     *
     *  @param   id - the institution id to search for
     *  @returns institution - matching institution for id or null
     */
    static async getInstitutionById(id: (ObjectId | string)): Promise<(IInstitution & {_id:ObjectId}) | null> {
        return await institutionModel.findById(id);
    }

    /**
     *  Updates the given institution corresponding to the id with the
     *  update_object.
     *
     *  @param institution_id
     *  @param update_object - includes fields to be updated
     */
    static async updateInstitution(institution_id: (ObjectId | string), update_object: UpdateInstitutionDTO) {
        await institutionModel.updateOne({_id: institution_id}, update_object);
    }
}
