import { IUser, userModel } from "../models/user";
import { AddUserDTO, UpdateUserDTO } from "../dtos/user.dto";
import { ObjectId, FilterQuery } from "mongoose";

/**
 *  @class UserService
 *
 *  Provides useful methods to access the database and modify users,
 *  which can be used by the route-controllers.
 */
export default class UserService {
  /**
   *  Gets all the users from the database.
   *
   *  @returns  allUsers
   */
  static async getAllUsers(): Promise<(IUser & { _id: ObjectId })[]> {
    return await userModel.find().exec();
  }

  /**
   * Search for all users by a keyword. Advanced search not implemented yet,
   * just matching first/last name at the moment
   * 
   * Does not return private info like email, password, creation date, etc.
   * 
   * Sorting methods (TODO):
   * default: ascending by name
   * "namedesc": descending by name
   *
   * @param keyword a keyword
   * @param sortBy sorting method
   * @returns array of users
   */
  static async searchUsers(keyword: string | null | undefined, sortBy: string | null | undefined) {
    //TODO implement actual fuzzy searching and other sorting methods
    let keywordFilter: FilterQuery<IUser> = {};
    let sort = {};
    sort[sortBy]=1;
    if (keyword && keyword.length > 0) {
      keywordFilter = { $or: [{ firstName: { $regex : "^" +  keyword, $options : 'i'}}, { lastName: { $regex : "^" +  keyword, $options : 'i'}}] };
    }

    return await userModel
      .find(keywordFilter, {
        _id: 1,
        email: 0,
        password: 0,
        isEmailVerified: 0,
        isAdministrator: 0,
        createdAt: 0,
        updatedAt: 0,
      })
      .sort(sort);
  }

  /**
   *  Search for a user with the given user id and return if found.
   *
   *  @param   user_id - the user id to search for
   *  @param   includePassword - if true, the returned user includes the password field
   *  @returns user - matched user to user_id or null
   */
  static async getUserById(
    user_id: ObjectId | string,
    includePassword: boolean = false
  ): Promise<(IUser & { _id: ObjectId }) | null> {
    return includePassword
      ? await userModel.findById(user_id).select("+password").exec()
      : await userModel.findById(user_id).exec();
  }

  /**
   *  Search for a user with the given email and return if found.
   *
   *  @param   email - the email field to search for
   *  @param   includePassword - if true, the returned user includes the password field
   *  @returns user - matched user to email or null
   */
  static async getUserByEmail(
    email: string,
    includePassword: boolean = false
  ): Promise<(IUser & { _id: ObjectId }) | null> {
    return includePassword
      ? await userModel.findOne({ email }).select("+password").exec()
      : await userModel.findOne({ email }).exec();
  }

  /**
   *  Gets all unauthorized users.
   *
   *  @returns users - unauthorized users or null
   */
  static async getUnauthUsers(): Promise<(IUser & { _id: ObjectId })[]> {
    return await userModel.find({ isAuthorized: false }).exec();
  }

  /**
   *  Adds given user to the database.
   *
   *  @param    user - the user to add to the db
   *  @returns  userAdded - the added user
   */
  static async addUser(user: AddUserDTO): Promise<IUser> {
    let userAdded: IUser | undefined = undefined;
    userAdded = await userModel.create(user);
    return userAdded;
  }

  /**
   *  Updates the given user corresponding to the user_id with the
   *  update_object.
   *
   *  @param user_id
   *  @param update_object - includes fields to be updated
   */
  static async updateUser(user_id: ObjectId | string, update_object: UpdateUserDTO) {
    await userModel.updateOne({ _id: user_id }, update_object);
  }

  static async unsetUserAvatar(user_id: ObjectId | string): Promise<string | null | undefined> {
    let old = await userModel.findByIdAndUpdate(user_id, { $unset: { avatarUrl: "" } });
    return old?.avatarUrl;
  }
}
