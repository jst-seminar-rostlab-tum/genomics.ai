import {Request} from "express";

export interface ExtRequest extends Request {
    isAuthenticated?: boolean;  // declare optional property "isAuthenticated"
    user_id?: number;           // declare optional property "user_id"
}