import { IToken, tokenModel } from "../models/token";
import { AddTokenDTO } from "../dtos/token.dto";
import { ObjectId } from "mongoose";

/**
 *  @class TokenService
 *
 *  Provides useful methods to access the database and modify tokens,
 *  which can be used by the route-controllers.
 */
export default class TokenService {
  /**
   *  Creates token for the given userId.
   *
   *  @param    token - with userId field
   *  @returns  token
   */
  static async addToken(token: AddTokenDTO): Promise<IToken> {
    return await tokenModel.create({ _userId: token._userId });
  }

  /**
   *  Search token by the given userId.
   *
   *  @param    userId
   *  @returns  token or null
   */
  static async getTokenByUserId(user_id: ObjectId): Promise<(IToken & { _id: ObjectId }) | null> {
    return await tokenModel.findOne({ _userId: user_id });
  }

  /**
   *  Search token by the given token.
   *
   *  @param    token
   *  @returns  token or null
   */
  static async getTokenByToken(token: string): Promise<(IToken & { _id: ObjectId }) | null> {
    return await tokenModel.findOne({ token });
  }
}
