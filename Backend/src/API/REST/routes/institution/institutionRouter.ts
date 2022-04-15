import express, {Router} from "express";
import {IInstitution, institutionModel} from "../../../../database/models/institution";
import { userModel } from "../../../../database/models/user";
import check_auth from "../../middleware/check_auth";

export default function institution_route() : Router {
    let router = express.Router();

    router
        .post("/institutions", check_auth(), async (req: any, res) => {
        
            const {name, country, profilePictureURL, backgroundPictureURL} = req.body;
            const admin_user_id = req.user_id;

            if (!(name && country && admin_user_id))
                return res.status(400).send("Missing parameters");

            if (await institutionModel.findOne({name}))
                return res.status(409).send("Institution with the given name already exists!");

            if (! await userModel.findOne({_id: admin_user_id}))
                return res.status(404).send("Admin that you are trying to assign does not exists!");   

            try{
                const institution = await institutionModel.create({
                    name: name,
                    country: country,
                    profilePictureURL,
                    backgroundPictureURL,
                    adminIds: [admin_user_id]
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
