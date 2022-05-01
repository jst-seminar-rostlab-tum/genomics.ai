import { ObjectId } from "mongoose";

/**
 *  Stores the raw data needed to create a token.
 */
export interface AddTokenDTO {
  _userId: ObjectId;
}
