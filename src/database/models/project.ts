import {Document, model, Schema} from "mongoose";

enum ProjectJobStatus {
    "UPLOAD_PENDING",
    "UPLOAD_COMPLETE",
    "PROCESSING_PENDING",
    "ABORTED",
    "DONE"
}

export interface IProject extends Document {
    owner: Schema.Types.ObjectId;

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
    owner: {type: Schema.Types.ObjectId, require: true},

    // file
    uploadId: {type: String, require: false},
    fileName: {type: String, require: false},
    location: {type: String, require: false},
    fileSize: {type: Schema.Types.Number, require: false, default: -1},
    uploadDate: {type: Schema.Types.Date, require: true},

    // project
    status: {type: String, require: true, enum: ProjectJobStatus},

    resultName: {type: String, require: false},
    resultSize: {type: Schema.Types.Number, require: false, default: -1}
});

export const projectModel = model<IProject>("Project", projectSchema);

