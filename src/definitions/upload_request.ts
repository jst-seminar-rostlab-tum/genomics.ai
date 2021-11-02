import {ObjectId} from "mongoose";
import {ExtRequest} from "./ext_request";

export interface UploadRequest extends ExtRequest {
    is_authenticated?: boolean;  // declare optional property "is_authenticated"
    user_id?: ObjectId;         // declare optional property "user_id"
    email?: string;              // declare optional property "email"
    uploadId: string;             // declare required property "uploadId"
}

