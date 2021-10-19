import {Document, Schema, model} from "mongoose";

export interface IUser extends Document {
    firstName: string;
    lastName: string;
    email: string;
    password: string;
    note: string;
}

const userSchema = new Schema<IUser>({
    firstName: {type: String, require: true},   // needed for contact-emails
    lastName: {type: String, default: ""},
    email: {type: String, unique: true, require: true},
    password: {type: String, require: true, select: false},
    token: {type: String}
});

export const userModel = model<IUser>("User", userSchema);
