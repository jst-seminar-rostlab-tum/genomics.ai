import { Document, Schema, model } from "mongoose";
import * as crypto from "crypto";

export interface IProjectUpdateToken extends Document {
  _projectId: Schema.Types.ObjectId;
  token: string;
}

const projectUpdateTokenSchema = new Schema<IProjectUpdateToken>({
  _projectId: { type: Schema.Types.ObjectId, required: true, ref: "Project" },
  token: { type: String, required: true, default: () => crypto.randomBytes(32).toString("hex") },
});

export const projectUpdateTokenModel = model<IProjectUpdateToken>(
  "ProjectUpdateToken",
  projectUpdateTokenSchema,
);