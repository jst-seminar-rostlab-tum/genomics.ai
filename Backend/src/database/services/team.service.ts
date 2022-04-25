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
     *  @param    project
     *  @returns  projectAdded - the added team
     */
    static async addTeam(project: AddTeamDTO): Promise<ITeam> {
        let projectAdded : (ITeam | undefined) = undefined;
        projectAdded = await teamModel.create(project);
        return projectAdded;
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
     *  @param   projectId
     *  @returns project - matched proejct to projectId or null
     */
    static async getTeamById(projectId: (ObjectId | string)):
      Promise<( ITeam & { _id: ObjectId } | null )> {
        return await teamModel.findById(projectId).exec();
    }

    /**
     *  Add the given userId to the invitation list of the given team.
     *
     *  @param   projectId
     *  @param   userId
     *  @returns updateDocument
     */
    static async addInvitationMemberId(projectId: (ObjectId | string), userId: (ObjectId | string)): Promise<any> {
        return await teamModel.updateOne(
            { _id: projectId },
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
}
