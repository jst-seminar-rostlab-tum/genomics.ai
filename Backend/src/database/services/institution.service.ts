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
