import {Document, model, Schema} from "mongoose";

export interface IProject extends Document {
    owner: string;
    location: string;
    fileName: string;
    uploadId: string;
}

const projectSchema = new Schema<IProject>({
    owner: {type: String, require: true},
    location: {type: String, require: true, default: ""},
    fileName: {type: String, require: true},
    uploadId: {type: String, require: true}
});

export const projectModel = model<IProject>("Project", projectSchema);
