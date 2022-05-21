import { Document, Schema, model } from "mongoose";
import * as crypto from "crypto";

export interface IPasswordResetToken extends Document {
  _userId: Schema.Types.ObjectId;
  token: string;
  createdAt: Date;
}

const passwordResetTokenSchema = new Schema<IPasswordResetToken>({
  createdAt: { type: Date, default: () => new Date(), expires: 10800 },
  _userId: { type: Schema.Types.ObjectId, required: true, ref: "User" },
  token: { type: String, required: true, default: () => crypto.randomBytes(32).toString("hex") },
});

export const passwordResetTokenModel = model<IPasswordResetToken>(
  "PasswordResetToken",
  passwordResetTokenSchema
);
