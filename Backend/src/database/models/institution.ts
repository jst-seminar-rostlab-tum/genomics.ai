import {Document, model, Schema} from "mongoose";

export interface IInstitution extends Document {
   name:string,
   country:string,
   profilePictureURL: string,
   backgroundPictureURL: string,
   adminIds: [Schema.Types.ObjectId],
   invitedMemberIds: [Schema.Types.ObjectId],
   projects: [Schema.Types.ObjectId],
}

const institutionSchema = new Schema<IInstitution>({
    name: {
        type: String,
        required: true
    },

    country: {
        type: String, 
        require: true
    },

    profilePictureURL: {
        type: String, 
        require: false
    },

    backgroundPictureURL: {
        type: String, 
        require: false
    },

    adminIds: [{
        type: Schema.Types.ObjectId, 
        ref: 'User',
        require: false
    }],

    invitedMemberIds: [{
        type: Schema.Types.ObjectId, 
        ref: 'User',
        require: false
    }],

    projects: [{
        type: Schema.Types.ObjectId, 
        ref: 'Project',
        require: false
    }],
    
}, {
    timestamps: true,
});

export const institutionModel = model<IInstitution>("Institution", institutionSchema);

