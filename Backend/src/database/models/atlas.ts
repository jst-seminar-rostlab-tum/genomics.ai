import {Document, model, Schema} from "mongoose";

export interface IAtlas extends Document {
   name: string,
   previewPictureURL: string,
   modalities: Array<string>,
   numberOfCells: number,
   species: Array<string>,
   compatibleModels: Array<Schema.Types.ObjectId>
}

const atlasSchema = new Schema<IAtlas>({
    name: {
        type: String,
        required: true
    },

    previewPictureURL: {
        type: String,
        required: false,
    },

    modalities: [{
        type: String,
        required: true
    }],

    numberOfCells: {
        type: Number,
        required: true
    },

    species: [{
        type: String,
        required: true
    }],

    compatibleModels: [{
        type: String,
        required: true
    }]

}, {
    timestamps: true,
});

export const atlasModel = model<IAtlas>("Atlas", atlasSchema);
