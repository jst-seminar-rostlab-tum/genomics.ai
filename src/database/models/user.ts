import {Document, Schema, model} from "mongoose";

export interface IUser extends Document {
    firstName: string;
    lastName: string;
    email: string;
    password: string;
    note: string;
    authorized: boolean;
    token: string
}

const userSchema = new Schema<IUser>({
    firstName: {type: String, require: true},   // needed for contact-emails
    lastName: {type: String, default: ""},
    email: {type: String, unique: true, require: true},
    password: {type: String, require: true, select: false},
    note: {type: String, require: false},
    authorized: {type: Boolean, required: true, default: false}
});

export const userModel = model<IUser>("User", userSchema);
