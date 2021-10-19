import {Request} from "express";
import {ObjectId} from "mongoose";

export interface ExtRequest extends Request {
    is_authenticated?: boolean;  // declare optional property "is_authenticated"
    user_id?: ObjectId ;         // declare optional property "user_id"
    email?: string;              // declare optional property "email"
}

