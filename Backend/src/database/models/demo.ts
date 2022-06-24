import { Document, model, Schema } from "mongoose";

export interface IDemo extends Document {
    name: string;
    model: string; 
    atlas: string;
    dataURL: string;
}

const demoSchema = new Schema<IDemo>({
    name: { type: String, required: true},
    model: { type: String, required: true},
    atlas: { type: String, required: true},
    dataURL: { type: String, required: true},
})

export const demoModel = model<IDemo>("Demo", demoSchema);


// a schema defines the structure of the object stored in it
// a model is the constructor for the schema
// a document is an instance of the model