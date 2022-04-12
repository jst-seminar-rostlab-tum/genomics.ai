import {Document, Schema, model} from "mongoose";

export interface IUser extends Document {
    firstName: string;
    lastName: string;
    email: string;
    password?: string;
    note: string;
    token: string;
    emailVerificationToken: string;
    isEmailVerified: boolean; // email has been verified
    isAuthorized: boolean; // user has been authorized by an administrator
    isAdministrator: boolean;
}

const userSchema = new Schema<IUser>({
    firstName: {type: String, required: true},   // needed for contact-emails
    lastName: {type: String, default: ""},
    email: {type: String, unique: true, required: true},
    password: {type: String, required: true, select: false},
    note: {type: String, required: false},
    isEmailVerified: {type: Boolean, required: true, default: false},
    isAdministrator: {type: Boolean, required: true, default: false}
});

export const userModel = model<IUser>("User", userSchema);
