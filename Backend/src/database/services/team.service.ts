import { ITeam, teamModel } from "../models/team";
import { IUser } from "../models/user";
import { AddTeamDTO, UpdateTeamDTO } from "../dtos/team.dto";
import { ObjectId } from "mongoose";

/**
 *  @class TeamService
 *
 *  Provides useful methods to access the database and modify teams,
 *  which can be used by the route-controllers.
 */
export default class TeamService {
  public static mergeAdminsMembers<
    T extends T2 | Array<T2>,
    T2 extends { memberIds: Array<any>; adminIds: Array<any> }
  >(team: T): T {
    if (Array.isArray(team)) {
      for (let t of team) {
        t.memberIds.push(...t.adminIds);
      }
    } else {
      team.memberIds.push(...team.adminIds);
    }
    return team;
  }

  /**
   *  Adds given team to the database.
   *
   *  @param    team
   *  @returns  teamAdded - the added team
   */
  static async addTeam(team: AddTeamDTO): Promise<ITeam> {
    let teamAdded: ITeam | undefined = undefined;
    teamAdded = await teamModel.create(team);
    return teamAdded;
  }

  static async updateTeam(updateTeam: UpdateTeamDTO): Promise<any> {
    const updateObj: any = {};
    if (updateTeam.description) updateObj.description = updateTeam.description;
    if (updateTeam.visibility) updateObj.visibility = updateTeam.visibility;
    const updateResult = await teamModel.updateOne({ _id: updateTeam.id }, { $set: updateObj });

    return updateResult;
  }

  /**
   *  Search for a team with the given title and return if found.
   *
   *  @param   title
   *  @returns team or null
   */
  static async getTeamByTitle(title: string): Promise<(ITeam & { _id: ObjectId }) | null> {
    return await teamModel.findOne({ title });
  }

  /**
   *  Search for a team with the given team id and return if found.
   *
   *  @param   teamId
   *  @returns team - matched team to teamId or null
   */
  static async getTeamById(teamId: ObjectId | string): Promise<(ITeam & { _id: ObjectId }) | null> {
    return await teamModel.findById(teamId).exec();
  }

  /**
   *  Return all the teams that has the member specified by userId.
   *
   *  @param   userId
   *  @returns array of teamId's and titles
   */
  static async getTeamsOfUser(userId: ObjectId): Promise<(ITeam & { _id: ObjectId })[]> {
    return await teamModel
      .find(
        {
          $or: [{ memberIds: userId }, { adminIds: userId }],
        },
        { title: 1 }
      )
      .exec();
  }

  static async getMembersOfTeam(
    team_id: ObjectId | string
  ): Promise<{ memberIds: Array<IUser>; adminIds: Array<IUser> } | null> {
    let team = await teamModel.findById(team_id).populate("memberIds").populate("adminIds");
    if (!team) return null;
    let { memberIds, adminIds } = team;
    return {
      memberIds: memberIds as any as Array<IUser>,
      adminIds: adminIds as any as Array<IUser>,
    };
    //Cast Array<ObjectId> to Array<IUser> as populate changes the type of elements of memberIds
  }

  /**
   *  Add the given userId to the invitation list of the given team.
   *
   *  @param   teamId
   *  @param   userId
   *  @returns updateDocument
   */
  static async addInvitationMemberId(
    teamId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await teamModel.updateOne({ _id: teamId }, { $addToSet: { invitedMemberIds: userId } });
  }

  /**
   *  Add the given userId to the admin list and removes he/she from the memberIds, of the given team.
   *
   *  @param   teamId
   *  @param   userId
   *  @returns updateDocument
   */
  static async addAdminToTeam(teamId: ObjectId | string, userId: ObjectId | string): Promise<any> {
    return await teamModel.updateOne(
      { _id: teamId },
      {
        $addToSet: { adminIds: userId },
        $pull: { memberIds: userId },
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
  static async isAdmin(
    userId: ObjectId | string,
    team: (ITeam & { _id: ObjectId }) | null
  ): Promise<boolean> {
    if (!team) return false; /* team does not exist */

    let isAdmin = false;
    var listAdmins = team.adminIds.map(String);
    var userIdStr = String(userId);
    if (listAdmins.includes(userIdStr)) isAdmin = true;
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
  static async isMember(
    userId: ObjectId | string,
    team: (ITeam & { _id: ObjectId }) | null
  ): Promise<boolean> {
    if (!team) return false; /* team does not exist */

    let isMember = false;
    var listMembers = team.memberIds.map(String);

    var userIdStr = String(userId);
    if (listMembers.includes(userIdStr)) isMember = true;
    return isMember;
  }

  /**
   *  Add the given userId into the given team.
   *
   *  @param   teamId
   *  @param   userId
   *  @returns updateDocument
   */
  static async joinMemberIntoTeam(
    teamId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await teamModel.updateOne(
      { _id: teamId },
      {
        $addToSet: { memberIds: userId },
        $pull: { invitedMemberIds: userId },
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
  static async setInstitutionOfTeam(
    teamId: ObjectId | string,
    institutionId: ObjectId | string
  ): Promise<any> {
    return await teamModel.updateOne({ _id: teamId }, { $set: { institutionId: institutionId } });
  }

  static async getTeams(queryParams: any): Promise<ITeam[] | null> {
    var filter: any,
      sortBy = {};

    queryParams.hasOwnProperty("keyword")
      ? (filter = { title: { $regex: "^" + queryParams.keyword, $options: "i" } })
      : (filter = {});
    queryParams.hasOwnProperty("visibility") ? (filter.visibility = queryParams.visibility) : null;

    if (queryParams.hasOwnProperty("sortBy")) {
      if (queryParams.sortBy == "title") sortBy["title"] = 1;
      else sortBy["updatedAt"] = -1;
    } else sortBy = {};

    return await teamModel.find(filter).populate("institutionId").sort(sortBy);
  }

  static async getUsersTeams(userId: ObjectId | string): Promise<ITeam[] | null> {
    return await teamModel.find({
      $or: [
        { memberIds: { $elemMatch: { $eq: userId } } },
        { adminIds: { $elemMatch: { $eq: userId } } },
      ],
    });
  }

  /**
   *  Add the given userId to the admin list and removes he/she from the memberIds, of the given team.
   *
   *  @param   teamId
   *  @param   institutionId
   *  @returns updateDocument
   */
  static async removeTeamFromInstitution(
    teamId: ObjectId | string,
    institutionId: ObjectId | string
  ): Promise<any> {
    return await teamModel.updateOne(
      { _id: teamId },
      {
        $unset: { institutionId: 1 },
      }
    );
  }

  /**
   *  Remove the given userId from the given team.
   *
   *  @param   teamId
   *  @param   userId
   *  @returns updateDocument
   */
  static async removeMemberFromTeam(
    teamId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await teamModel.updateOne(
      { _id: teamId },
      {
        $pull: { memberIds: userId, adminIds: userId },
      }
    );
  }

  /**
   *  Remove the given userId from the admins list of the given team.
   *
   *  @param   teamId
   *  @param   userId
   *  @returns updateDocument
   */
  static async demoteAdminFromTeam(
    teamId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await teamModel.updateOne(
      { _id: teamId },
      {
        $pull: { adminIds: userId },
        $push: { memberIds: userId },
      }
    );
  }
}
