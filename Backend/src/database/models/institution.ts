import { Document, model, Schema } from "mongoose";

export enum visibilityInstitutionStatus {
  "PRIVATE",
  "PUBLIC",
}

export interface IInstitution extends Document {
  name: string;
  country: string;
  profilePictureURL: string;
  backgroundPictureURL: string;
  adminIds: Array<Schema.Types.ObjectId>;
  memberIds: Array<Schema.Types.ObjectId>;
  invitedMemberIds: Array<Schema.Types.ObjectId>;
  visibility: string;
}

const institutionSchema = new Schema<IInstitution>(
  {
    name: {
      type: String,
      required: true,
    },

    country: {
      type: String,
      require: true,
    },

    profilePictureURL: {
      type: String,
      require: false,
    },

    backgroundPictureURL: {
      type: String,
      require: false,
    },

    adminIds: [
      {
        type: Schema.Types.ObjectId,
        ref: "User",
        require: true,
      },
    ],

    memberIds: [
      {
        type: Schema.Types.ObjectId,
        ref: "User",
        require: false,
      },
    ],

    invitedMemberIds: [
      {
        type: Schema.Types.ObjectId,
        ref: "User",
        require: false,
      },
    ],

    visibility: {
      type: String,
      enum: visibilityInstitutionStatus,
      require: true,
    },
  },
  {
    timestamps: true,
  }
);

export const institutionModel = model<IInstitution>("Institution", institutionSchema);
