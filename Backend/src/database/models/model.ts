import { Document, model, Schema } from "mongoose";

export interface IModel extends Document {
  name: string;
  description: string;
  requirements: Array<string>;
}

const modelSchema = new Schema<IModel>(
  {
    name: {
      type: String,
      required: true,
    },

    description: {
      type: String,
      required: false,
      default: "",
    },

    requirements: [
      {
        type: Schema.Types.ObjectId,
        required: false,
        default: "",
      },
    ],
  },
  {
    timestamps: true,
  }
);

export const modelModel = model<IModel>("Model", modelSchema);
