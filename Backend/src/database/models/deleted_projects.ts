import { Document, model, Schema } from "mongoose";
import { IProject, projectModel } from "./project";

export const DELETED_PROJECT_LIFETIME_DAYS = Number.parseFloat(process.env.PROJECT_RECYCLE_BIN_LIFETIME_DAYS);

export interface IDeletedProject extends IProject {
  deletedAt: Date;
}

const projectSchema = projectModel.schema.obj;
const extendedSchema = Object.assign({}, projectSchema, {
  deletedAt: { type: Schema.Types.Date, require: true },
});

const deletedProjectSchema = new Schema<IDeletedProject>(extendedSchema);

export const deletedProjectModel = model<IDeletedProject>("DeletedProject", deletedProjectSchema);
