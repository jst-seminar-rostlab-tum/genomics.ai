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
    static async addTeam(team: AddTeamDTO): Promise<ITeam> {
        let teamAdded : (ITeam | undefined) = undefined;
        teamAdded = await teamModel.create(team);
        return teamAdded;
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
     *  Return all the teams that has the member specified by userId.
     *
     *  @param   userId
     *  @returns array of teamId's and titles
     */
    static async getTeamsOfUser(userId: (ObjectId)):
      Promise<( ITeam & { _id: ObjectId } )[] > {
        return await teamModel.find( {
          $or: [ { memberIds: userId }, { adminIds: userId } ] },
          { title: 1 } ).exec();
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
     *  Add the given projectId to the project list of the given team.
     *
     *  @param   teamId
     *  @param   projectId
     *  @returns updateDocument
     */
    static async addProject(teamId: (ObjectId | string), projectId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { $addToSet: { projects: projectId} }
        );
    }

    /**
     *  Add the given userId to the admin list and removes he/she from the memberIds, of the given team.
     *
     *  @param   teamId
     *  @param   userId
     *  @returns updateDocument
     */
    static async addAdminToTeam(teamId: (ObjectId | string), userId: (ObjectId | string)):
      Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            {
                $addToSet: { adminIds: userId },
                $pull: { memberIds: userId }
            }
        );
    }

    /**
     *  Returns true if the given user is an admin of the given team.
     *  The given team should exist, otherwise the method returns false.
     *
     *  @param  userId
     *  @param  team
     *  @return isAdmin
     */
    static async isAdmin(userId: (ObjectId | string), team: (ITeam & {_id: ObjectId}) | null):
      Promise<boolean> {
        if (!team)
          return false; /* team does not exist */

        let isAdmin = false;
        var listAdmins = team.adminIds.map(String);
        var userIdStr= String(userId);
        if (listAdmins.includes(userIdStr))
            isAdmin = true;
        return isAdmin;
    }

    /**
     *  Returns true if the given user is a member of the given team.
     *  The given team should exist, otherwise the method returns false.
     *
     *  @param  userId
     *  @param  team
     *  @return isMember
     */
    static async isMember(userId: (ObjectId | string), team: (ITeam & {_id: ObjectId}) | null):
      Promise<boolean> {
        if (!team)
          return false; /* team does not exist */

        let isMember = false;
        var listMembers = team.memberIds.map(String);

        var userIdStr = String(userId);
        if (listMembers.includes(userIdStr))
            isMember = true;
        return isMember;
    }

    /**
     *  Add the given userId into the given project.
     *
     *  @param   teamId
     *  @param   userId
     *  @returns updateDocument
     */
    static async addNewMemberIntoTeam(teamId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: teamId },
            { $addToSet: { memberIds: userId} }
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
}
