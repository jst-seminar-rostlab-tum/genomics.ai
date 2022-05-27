import { ObjectId } from "mongoose";

export function query_path(projectid: ObjectId | string): string {
  return `projects/${projectid}/query.h5ad`;
}

export function result_path(projectid: ObjectId | string): string {
  return `results/${projectid}/query.csv`;
}
export function result_model_path(projectid: ObjectId | string): string {
  return `results/${projectid}/model.pt`;
}
