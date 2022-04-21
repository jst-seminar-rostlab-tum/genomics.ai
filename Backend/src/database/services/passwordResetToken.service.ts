import {IPasswordResetToken, passwordResetTokenModel} from "../models/password_reset_token";
import {AddPasswordResetTokenDTO} from "../dtos/password_reset_token.dto";
import {ObjectId} from "mongoose";

/**
 *  @class PasswordResetTokenService
 *
 *  Provides useful methods to access the database and modify
 *  password-reset-tokens, which can be used by the route-controllers.
 */
export default class PasswordResetTokenService {
    /**
     *  Creates password-reset-token for the given userId.
     *
     *  @param    password-reset-token - with userId field
     *  @returns  password-reset-token
     */
    static async addToken(token: AddPasswordResetTokenDTO): Promise<IPasswordResetToken> {
        return await passwordResetTokenModel.create({ _userId: token._userId });
    }

    /**
     *  Search password-reset-token by the given token.
     *
     *  @param    token
     *  @returns  password-reset-token or null
     */
    static async getTokenByToken(token: string):
      Promise<( IPasswordResetToken & { _id: ObjectId } | null )> {
        return await passwordResetTokenModel.findOne({token});
    }
}
