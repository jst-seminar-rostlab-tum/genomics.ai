import {Document, model, Schema} from "mongoose";

export enum ProjectStatus {
    UPLOAD_PENDING="UPLOAD_PENDING",
    UPLOAD_COMPLETE="UPLOAD_COMPLETE",
    PROCESSING_PENDING="PROCESSING_PENDING",
    PROCESSING_FAILED="PROCESSING_FAILED",
    ABORTED="ABORTED",
    DONE="DONE",
}

export interface IProject extends Document {
    owner: Schema.Types.ObjectId;
    teamId: Schema.Types.ObjectId;
    name: string;

    modelId: string;
    atlasId: string;

    // file
    uploadId: string;
    location: string;
    fileName: string;
    fileSize: number; // (of bytes)
    uploadDate: Date;

    // project
    status: string;
    resultName: string;
    resultSize: number;
}

const projectSchema = new Schema<IProject>({
    owner: {
        type: Schema.Types.ObjectId, 
        ref: 'User',
        require: true
    },
    teamId: {
        type: Schema.Types.ObjectId,
        ref: 'Team',
        required: false
    },
    name: {type: String, require: true},

    modelId: {type: String, require: true},
    atlasId: {type: String, require: true},

    // file
    uploadId: {type: String, require: false},
    fileName: {type: String, require: false},
    location: {type: String, require: false},
    fileSize: {type: Schema.Types.Number, require: false, default: -1},
    uploadDate: {type: Schema.Types.Date, require: true},

    // project
    status: {type: String, require: true, enum: ProjectStatus},

    resultName: {type: String, require: false},
    resultSize: {type: Schema.Types.Number, require: false, default: -1}
});

export const projectModel = model<IProject>("Project", projectSchema);

