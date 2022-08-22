import express, { Router } from "express";
import { AddUserDTO } from "../../../database/dtos/user.dto";
import { IUser } from "../../../database/models/user";
import UserService from "../../../database/services/user.service";
import jwt from "jsonwebtoken";

const JWT_SECRET = process.env.JWT_SECRET || "";

/**
 * Temporary authentication route for the non-login version.
 * To be used as a one-time access for creating a project and processing it. 
 */
export default function get_temp_auth(): Router {
    let router = express.Router();
    
    router.get("/temp_auth", async (req: any, res) => {
        let tempUser: IUser | undefined = undefined;
        // Create temporary user that will be deleted after 
        // upload/download process is finished.  
        try {
            // TODO: delete user afterwards  
            let userToAdd: AddUserDTO = {
                firstName: "-",
                lastName: "-",
                email: "-",
                password: "-",
                note: "temporary_user"
            };
            // add temporary user
            tempUser = await UserService.addUser(userToAdd);

            // update the email with the unique id value
            // in order to bypass the uniqueness requirement of the email
            UserService.updateUser(tempUser._id, { email: tempUser._id })

            // create JWT token for temporary access, use the temporary id for both the email and id.
            const token = jwt.sign({ id: tempUser._id, email: tempUser._id }, JWT_SECRET, {
            expiresIn: "20h",
            });

            // return successful temp access response
            return res.status(200).json({
                msg: "Temporary access granted",
                user: tempUser.toObject(),
                jwt: token,
            });

        } catch (err) {
            console.error("Error registering user!");
            console.error(JSON.stringify(err));
            return res.status(500).send("Unknown error in database");
        }
    });
    
    return router;
}