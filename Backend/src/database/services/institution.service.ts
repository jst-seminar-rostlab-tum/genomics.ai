import { IInstitution, institutionModel } from "../models/institution";
import { IUser } from "../models/user";
import {ITeam, teamModel} from "../models/team";
import ProjectService from "./project.service";
import { AddInstitutionDTO, UpdateInstitutionDTO } from "../dtos/institution.dto";
import { ObjectId } from "mongoose";
import { IProject, projectModel } from "../models/project";
import TeamService from "./team.service";

/**
 *  @class InstitutionService
 *
 *  Provides useful methods to access the database and modify institutions,
 *  which can be used by the route-controllers.
 */
export default class InstitutionService {
  public static mergeAdminsMembers<
    T extends T2 | Array<T2>,
    T2 extends { memberIds: Array<any>; adminIds: Array<any> }
  >(institution: T): T {
    if (Array.isArray(institution)) {
      for (let i of institution) {
        i.memberIds.push(...i.adminIds);
      }
    } else {
      institution.memberIds.push(...institution.adminIds);
    }
    return institution;
  }
  /**
   *  Adds given institution to the database.
   *
   *  @param    institution
   *  @returns  institutionAdded - the added institution
   */
  static async addInstitution(institution: AddInstitutionDTO): Promise<IInstitution> {
    let institutionAdded: IInstitution | undefined = undefined;
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
  static async inviteToInstitution(
    institutionId: ObjectId,
    userId: ObjectId
  ): Promise<IInstitution | undefined> {
    let updatedInstitution: IInstitution | undefined = undefined;

    const institution = await institutionModel.findOne({ _id: institutionId });

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
  static async makeUserAnAdminOfInstitution(
    institutionId: ObjectId,
    userId: ObjectId
  ): Promise<IInstitution | undefined> {
    let updatedInstitution: IInstitution | undefined = undefined;

    const institution = await institutionModel.findOne({ _id: institutionId });

    if (institution) {
      institution.adminIds = [...institution.adminIds, userId];
      institution.memberIds = institution.memberIds.filter((id) => String(id) !== String(userId));
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
  static async makeUserMemberOfInstitution(
    institutionId: ObjectId,
    userId: ObjectId
  ): Promise<IInstitution | undefined> {
    let updatedInstitution: IInstitution | undefined = undefined;

    const institution = await institutionModel.findOne({ _id: institutionId });

    if (institution) {
      institution.memberIds = [...institution.memberIds, userId];
      institution.invitedMemberIds = institution.invitedMemberIds.filter(
        (id) => String(id) !== String(userId)
      );
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
  static async findMemeberOrInvitedById(
    user_id: ObjectId | string,
    institution_id: ObjectId | string
  ): Promise<(IInstitution & { _id: any }) | undefined> {
    const result = await institutionModel.findOne({
      _id: institution_id,
      $or: [
        {
          invitedMemberIds: { $elemMatch: { $eq: user_id } },
        },
        {
          memberIds: { $elemMatch: { $eq: user_id } },
        },
      ],
    });
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
  static async findMemeberById(
    user_id: ObjectId | string,
    institution_id: ObjectId | string
  ): Promise<(IInstitution & { _id: any }) | undefined> {
    const result = await institutionModel.findOne({
      _id: institution_id,
      $or: [
        {
          memberIds: { $elemMatch: { $eq: user_id } },
        },
      ],
    });

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
  static async getInstitutionByName(
    name: string
  ): Promise<(IInstitution & { _id: ObjectId }) | null> {
    return await institutionModel.findOne({ name });
  }

  /**
   *  Search for an institution with the given id and return if found.
   *
   *  @param   id - the institution id to search for
   *  @returns institution - matching institution for id or null
   */
  static async getInstitutionById(
    institutionId: ObjectId | string
  ): Promise<(IInstitution & { _id: ObjectId }) | null> {
    return await institutionModel.findById(institutionId).exec();
  }

  /**
   *  Add the given teamId into the institution.
   *
   *  @param   teamId
   *  @param   institutionId
   *  @returns updateDocument
   */
  static async addNewMemberIntoTeam(
    teamId: ObjectId | string,
    institutionId: ObjectId | string
  ): Promise<any> {
    return await institutionModel.updateOne(
      { _id: institutionId },
      { $addToSet: { memberIds: teamId } }
    );
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
    return await institutionModel.updateOne(
      { _id: institutionId },
      {
        $pull: { memberIds: teamId },
      }
    );
  }

  /**
   *  Returns true if the given user is an admin of the given institution.
   *  The given institution should exist, otherwise the method returns false.
   *
   *  @param  userId
   *  @param  institution
   *  @return isAdmin
   */
  static async isAdmin(
    userId: ObjectId | string,
    institution: (IInstitution & { _id: ObjectId }) | null
  ): Promise<boolean> {
    if (!institution) return false; /* institution does not exist */

    let isAdmin = false;
    var listAdmins = institution.adminIds.map(String);
    var userIdStr = String(userId);
    if (listAdmins.includes(userIdStr)) isAdmin = true;
    return isAdmin;
  }

  /**
   *  Returns true if the given user is a member of the given institution.
   *  The given institution should exist, otherwise the method returns false.
   *
   *  @param  userId
   *  @param  institution
   *  @return isMember
   */
  static async isMember(
    userId: ObjectId | string,
    institution: (IInstitution & { _id: ObjectId }) | null
  ): Promise<boolean> {
    if (!institution) return false; /* institution does not exist */

    let isMember = false;
    var listMembers = institution.memberIds.map(String);

    var userIdStr = String(userId);
    if (listMembers.includes(userIdStr)) isMember = true;
    return isMember;
  }

  /*
   *  Updates the given institution corresponding to the id with the
   *  update_object.
   *
   *  @param institution_id
   *  @param update_object - includes fields to be updated
   */
  static async updateInstitution(
    institution_id: ObjectId | string,
    update_object: UpdateInstitutionDTO
  ) {
    await institutionModel.updateOne({ _id: institution_id }, update_object);
  }

  static async unsetProfilePicture(
    institution_id: ObjectId | string
  ): Promise<string | null | undefined> {
    let old = await institutionModel.findByIdAndUpdate(institution_id, {
      $unset: { profilePictureURL: "" },
    });
    return old?.profilePictureURL;
  }

  static async unsetBackgroundPicture(
    institution_id: ObjectId | string
  ): Promise<string | null | undefined> {
    let old = await institutionModel.findByIdAndUpdate(institution_id, {
      $unset: { backgroundPictureURL: "" },
    });
    return old?.backgroundPictureURL;
  }

  static async filterInstitutions(query: any): Promise<any | null> {
    var keyword: object, sortBy =  {};

    query.hasOwnProperty("keyword")
      ? (keyword = { name: { $regex: "^" + query.keyword, $options: "i" } })
      : (keyword = {});

    if (query.hasOwnProperty("sortBy")) {
      if(query.sortBy == "name")
        sortBy["name"] = 1;
      else
        sortBy["updatedAt"] = -1;
    } else sortBy = {};

    const institutions : any =  await institutionModel.find(keyword)
        .sort(sortBy).lean();

    var populatedInstitutions : Array<any> = new Array<any>();

    for( let institution of institutions){
      let institutionsTeams = await teamModel.find({ institutionId: institution._id } );
      institution = {...institution , teams: institutionsTeams}
      populatedInstitutions.push(institution);
    }

    return populatedInstitutions;
  }

  static async getMembersOfInstitution(
    institution_id: ObjectId | string
  ): Promise<{memberIds: Array<IUser>, adminIds: Array<IUser>} | null> {
    let institution =  await institutionModel.findById(institution_id).populate("memberIds").populate("adminIds");
    if(!institution) return null;
    let {memberIds, adminIds} = institution;
    return {memberIds: memberIds as any as Array<IUser>, adminIds: adminIds as any as Array<IUser>};
  }

  static async getProjectsOfInstitution(institution_id: ObjectId | string) {
    const institution = await institutionModel.findById(institution_id);

    if (!institution) return null;

    let projectsOfMembers: IProject[] = [];
    if (institution?.memberIds?.length > 0) {
      const resp = await ProjectService.getProjectsOfUsers(institution.memberIds);
      if (resp) {
        projectsOfMembers = resp;
      }
    }

    const teams = await teamModel.find({ institution_id: institution_id });
    const teamIds = teams.map((team) => team._id);
    const projectsOfTeams = await ProjectService.getProjectsOfTeams(teamIds);

    const projects = projectsOfMembers;
    for (const projOfTeam of projectsOfTeams) {
      const proj = projOfTeam as any;
      if (projects.find((p) => String(p._id) == String(proj._id))) {
        continue;
      }
      projects.push(proj);
    }

    return projects;
  }

  static async getTeamsOfInstitution(institutionId: ObjectId | any): Promise<ITeam[] | null> {
    return await teamModel.find({ institutionId: institutionId } );
  }


  static async getUsersInstitutions(userId: ObjectId | string): Promise<IInstitution[] | null> {
    return await institutionModel.find({
      $or: [
        { memberIds: { $elemMatch: { $eq: userId } } },
        { adminIds: { $elemMatch: { $eq: userId } } },
      ],
    });
  }

  /**
   *  Remove the given userId from the given institution.
   *
   *  @param   institutionId
   *  @param   userId
   *  @returns updateDocument
   */
  static async removeMemberFromTeam(
    institutionId: ObjectId | string,
    userId: ObjectId | string
  ): Promise<any> {
    return await institutionModel.updateOne(
      { _id: institutionId },
      {
        $pull: { memberIds: userId, adminIds: userId },
      }
    );
  }
}
