import {Document, Schema, model} from "mongoose";

export interface IUser extends Document {
    firstName: string;
    lastName: string;
    email: string;
    password: string;
    note: any;
    isEmailVerified: boolean; // Only for e-mail, not for academic affiliation
    isAccountApproved: boolean;
}

const userSchema = new Schema<IUser>({
    firstName: {type: String, require: true},   // needed for contact-emails
    lastName: {type: String, default: ""},
    email: {type: String, unique: true, require: true},
    password: {type: String, require: true, select: false},
    note: {type: Object, require: false},
    isEmailVerified: {type: Boolean, default: false},
    isAccountApproved: {type: Boolean, default: false}
});

export const userModel = model<IUser>("User", userSchema);
