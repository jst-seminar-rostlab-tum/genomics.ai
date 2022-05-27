import  { Document, Schema, model } from "mongoose";
import * as crypto from "crypto";

export enum JoinTarget {
  INSTITUTION = "INSTITUTION",
  TEAM = "TEAM",
}

export interface IJoinToken extends Document {
  target: JoinTarget,
  institutionId: Schema.Types.ObjectId;
  teamId: Schema.Types.ObjectId;
  userId: Schema.Types.ObjectId;
  token: string;
}

const joinTokenSchema = new Schema<IJoinToken>({
  target: {
    type: String,
    require: true,
  },
  userId: {
    type:Schema.Types.ObjectId,
    ref: "User",
    required: true
  },
  institutionId: {
    type: Schema.Types.ObjectId,
    ref: "Institution",
    required: institutionIdRequired,
  },
  teamId: {
    type: Schema.Types.ObjectId,
    ref: "Team",
    required: teamIdRequired,
  },
  token: {
    type: String,
    required: true,
    default: () => crypto.randomBytes(32).toString("hex")
  },
});

export const joinTokenModel = model<IJoinToken>(
  "JoinToken",
  joinTokenSchema
)

function institutionIdRequired() {
  return this.target==JoinTarget.INSTITUTION;
}
function teamIdRequired() {
  return this.target==JoinTarget.TEAM;
}
