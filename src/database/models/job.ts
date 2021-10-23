import {Document, Schema, model} from "mongoose";

export interface IJob extends Document {
    bucketId: number,
    isFinished: boolean,
}

const jobSchema = new Schema<IJob>({
    bucketId: {type: Number, require: true},
    isFinished: {type: Boolean, default: false},
});

export const jobModel = model<IJob>("Job", jobSchema);
