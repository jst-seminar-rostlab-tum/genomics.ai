import express, {Router} from "express";

import {Schema} from "mongoose";

import InstitutionService from "../../../../database/services/institution.service";
import {AddInstitutionDTO} from "../../../../database/dtos/institution.dto";
import UserService from "../../../../database/services/user.service";

import check_auth from "../../middleware/check_auth";

const create_institution = () : Router => {
    let router = express.Router();

    router
        .post("/institutions", check_auth(), async (req: any, res) => {

            const {name, country, profilePictureURL, backgroundPictureURL} = req.body;
            const admin_user_id = req.user_id;

            if (!(name && country && admin_user_id))
                return res.status(400).send("Missing parameters");

            const institution = await InstitutionService.getInstitutionByName(name);
            if (institution)
                return res.status(409).send("Institution with the given name already exists!");

            const user = await UserService.getUserById(admin_user_id);
            if (!user)
                return res.status(404).send("Admin that you are trying to assign does not exists!");

            try{
                const institutionToAdd: AddInstitutionDTO = {
                    name,
                    country,
                    profilePictureURL,
                    backgroundPictureURL,
                    adminIds: [admin_user_id]
                };
                const institution = await InstitutionService.addInstitution(institutionToAdd);

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

const invite_to_institution = () : Router => {
    let router = express.Router();

    router
        .put("/institutions/:id/invite", check_auth(), async (req: any, res) => {

            const { userId }: {userId: Schema.Types.ObjectId} = req.body;
            const institutionId_to_modify = req.params.id;

            try {

                if (!(userId))
                return res.status(400).send("Missing parameter");

                if (! await UserService.getUserById(userId))
                    return res.status(404).send("User that you are trying to invite does not exists!");      

                if (await InstitutionService.findMemeberById(userId,institutionId_to_modify))
                    return res.status(404).send("User that you are trying to invite to this institution already is an invited member or is a member!");      

                const updatedInstitution = await  InstitutionService.inviteToInstitution(institutionId_to_modify, userId)

                if(updatedInstitution) {
                    res.json(updatedInstitution);
                } else {
                    return res.status(409).send("Could not invite person to institution!");
                }

            } catch (error) {
                return res.status(500).send("Something went wrong: "+ error)
            }

        })

    return router;
}

const test_institution = () : Router => {
    let router = express.Router();

    router
        .get("/institutions/test", check_auth(), async (req: any, res) => {
            res.json("test router")
        })

    return router;
}

const get_institution = () : Router => {
    let router = express.Router();
    router
        .get("/institution/:id", check_auth(), async (req: any, res) =>{
            const institutionId = req.params.id;
            try {
                const institution = await InstitutionService.getInstitutionById(institutionId);
                return res.status(200).json(institution);
            } catch (err) {
                console.error(JSON.stringify(err));
                return res.status(404).send(`Institution ${institutionId} not found`);
            }
        })
    return router;
}

const get_institutions = () : Router => {
    let router = express.Router();
    router
        .get("/institutions", check_auth(), async( req: any, res) => {
            const query = {...req.query };
            try{
                const institutions = await InstitutionService.filterInstitutions(query);
                return res.status(200).json(institutions);
            } catch (err){
                console.error(JSON.stringify(err));
                return res.status(404).send(`No institutions found`);
            }
        })
    return router;
}

export { create_institution, test_institution, invite_to_institution, get_institutions, get_institution }
