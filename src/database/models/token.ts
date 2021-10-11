import {Document, Schema, model} from "mongoose";
import * as crypto from "crypto";

interface IToken extends Document {
    _userId: Schema.Types.ObjectId;
    token: string;
    createdAt: Date;
}

const tokenSchema = new Schema<IToken>({
    _userId: { type: Schema.Types.ObjectId, required: true, ref: 'User'},
    token: { type: String, default: () => crypto.randomBytes(16).toString("hex") },
    createdAt: { type: Date, expires: 86400, default: () => new Date() }
});

export const tokenModel = model<IToken>("Token", tokenSchema);