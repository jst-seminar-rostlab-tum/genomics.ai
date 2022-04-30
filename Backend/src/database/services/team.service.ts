import {ITeam, teamModel} from "../models/team";
import {AddTeamDTO} from "../dtos/team.dto";
import {ObjectId} from "mongoose";

/**
 *  @class TeamService
 *
 *  Provides useful methods to access the database and modify projects,
 *  which can be used by the route-controllers.
 */
export default class TeamService {
    /**
     *  Adds given team to the database.
     *
     *  @param    team
     *  @returns  projectAdded - the added team
     */
    static async addTeam(project: AddTeamDTO): Promise<ITeam> {
        let team : (ITeam | undefined) = undefined;
        team = await teamModel.create(team);
        return team;
    }

    /**
     *  Search for a team with the given title and return if found.
     *
     *  @param   title
     *  @returns project or null
     */
    static async getTeamByTitle(title: string):
      Promise<( ITeam & { _id: ObjectId } | null )> {
        return await teamModel.findOne({title});
    }

    /**
     *  Search for a team with the given team id and return if found.
     *
     *  @param   teamId
     *  @returns project - matched proejct to projectId or null
     */
    static async getTeamById(teamId: (ObjectId | string)):
      Promise<( ITeam & { _id: ObjectId } | null )> {
        return await teamModel.findById(teamId).exec();
    }

    /**
     *  Add the given userId to the invitation list of the given team.
     *
     *  @param   teamId
     *  @param   userId
     *  @returns updateDocument
     */
    static async addInvitationMemberId(teamId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { $addToSet: { invitedMemberIds: userId} }
        );
    }

    /**
     *  Add the given userId to the admin list and removes he/she from the memberIds, of the given team.
     *
     *  @param   teamId
     *  @param   userId
     *  @returns updateDocument
     */
     static async addAdminToTeam(teamId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { 
                $addToSet: { adminIds: userId },
                $pull: { memberIds: userId }
            }
        );
    }

    /**
     *  Add the given userId into the given project.
     *
     *  @param   teamId
     *  @param   userId
     *  @returns updateDocument
     */
    static async joinMemberIntoTeam(teamId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { 
                $addToSet: { memberIds: userId},
                $pull: { invitedMemberIds: userId }
            }
        );
    }

    /**
     *  Remove the given userId into the given project.
     *
     *  @param   teamId
     *  @param   userId
     *  @returns updateDocument
     */
     static async removeMemberFromTeam(teamId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { 
                $pull: { memberIds: userId, adminIds: userId}
            }
        );
    }

    /**
     *  Add the given teamId into the institution.
     *
     *  @param   teamId
     *  @param   institutionId
     *  @returns updateDocument
     */
    static async setInstitutionOfTeam(teamId: (ObjectId | string), institutionId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { $set : { institutionId: institutionId }}
        );
    }

    static async getTeams(queryParams: any ):
    Promise<ITeam[] | null >{

        var filter  : any,
            sortBy  : any;

        queryParams.hasOwnProperty('keyword') ?  filter = { name : queryParams.keyword } : filter = {};
        queryParams.hasOwnProperty('visibility') ? filter.visibility = queryParams.visibility : null;
        
        if(queryParams.hasOwnProperty('sortBy')){
            let sortProperty = queryParams.sortBy;
            sortBy = { sortProperty : 1 }
        } else
            sortBy = {};

        return await teamModel.find(filter).sort(sortBy);
    }

    /**
     *  Add the given userId to the admin list and removes he/she from the memberIds, of the given team.
     *
     *  @param   teamId
     *  @param   institutionId
     *  @returns updateDocument
     */
     static async removeTeamFromInstitution(teamId: (ObjectId | string), institutionId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            {
                $unset: { institutionId: 1 }
            }
        );
    }
}
