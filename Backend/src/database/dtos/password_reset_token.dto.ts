import { ObjectId } from "mongoose";

/**
 *  Stores the raw data needed to create a password reset token.
 */
export interface AddPasswordResetTokenDTO {
  _userId: ObjectId;
}
