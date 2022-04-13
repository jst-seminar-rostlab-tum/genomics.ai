import express, {Router} from "express";
import {IInstitution, institutionModel} from "../../../../database/models/institution";
import { userModel } from "../../../../database/models/user";

export default function institution_route() : Router {
    let router = express.Router();

    router
        .post("/institutions", async (req: any, res) => {
            
            const {name, country, profilePictureURL, backgroundPictureURL, adminId, invitedMemberIds, projects} = req.body;

            if (!(name && country && adminId))
                return res.status(400).send("Missing parameters");

            if (await institutionModel.findOne({name}))
                return res.status(409).send("Institution with the given name already exists!");

            if (! await userModel.findOne({_id: adminId}))
                return res.status(404).send("Admin that you are trying to assign does not exists!");   

            let institution : (IInstitution | undefined) = undefined;
            try{
                institution = await institutionModel.create({
                    name: name,
                    country: country,
                    profilePictureURL,
                    backgroundPictureURL,
                    invitedMemberIds,
                    projects,
                    adminIds: [adminId]
                });

                return res.status(201).json(institution);
            }catch(err){
                console.error("Error registering institution!");
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to create institution. (DB-error)");
            }
        })

    return router;
}
