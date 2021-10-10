import {Request} from "express";

export interface ExtRegisterRequest extends Request {
    body: {
        firstName: string;
        lastName: string;
        email: string;
        password: string;
    }
}