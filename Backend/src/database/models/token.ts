import { Document, Schema, model } from "mongoose";
import * as crypto from "crypto";

export interface IToken extends Document {
  _userId: Schema.Types.ObjectId;
  token: string;
  createdAt: Date;
}

const tokenSchema = new Schema<IToken>({
  createdAt: { type: Date, default: () => new Date() },
  _userId: { type: Schema.Types.ObjectId, required: true, ref: "User" },
  token: { type: String, required: true, default: () => crypto.randomBytes(16).toString("hex") },
});

tokenSchema.index({ createdAt: 1 }, { expireAfterSeconds: 86400 });

export const tokenModel = model<IToken>("Token", tokenSchema);
