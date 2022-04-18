import {Document, model, Schema} from "mongoose";

export enum visibilityStatus {
    "PRIVATE",
    "PUBLIC",
    "BY_INSTITUTION",
}
export interface IProject extends Document {
   title:string,
   description:string,
   adminIds: [Schema.Types.ObjectId],
   invitedMemberIds: [Schema.Types.ObjectId],
   memberIds: [Schema.Types.ObjectId],
   visibility: string,
   institutionId: Schema.Types.ObjectId,
}

const projectSchema = new Schema<IProject>({
    title: {
        type: String,
        required: true
    },

    description: {
        type: String, 
        require: true
    },

    adminIds: [{
        type: Schema.Types.ObjectId, 
        ref: 'User',
        require: true
    }],

    invitedMemberIds: [{
        type: Schema.Types.ObjectId, 
        ref: 'User',
        require: false
    }],

    memberIds: [{
        type: Schema.Types.ObjectId, 
        ref: 'User',
        require: false
    }],

    visibility: {
        type: String,
        enum: visibilityStatus,
        required: true
    },

    institutionId: {
        type: Schema.Types.ObjectId, 
        ref: 'Institution',
        require: false
    }
    
}, {
    timestamps: true,
});

export const projectModel = model<IProject>("Project", projectSchema);

