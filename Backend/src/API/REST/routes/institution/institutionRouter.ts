import express, {Router} from "express";
import { institutionModel} from "../../../../database/models/institution";
import { userModel } from "../../../../database/models/user";
import check_auth from "../../middleware/check_auth";

const create_institution = () : Router => {
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

const test_institution = () : Router => {
    let router = express.Router();

    router
        .get("/institutions/test", check_auth(), async (req: any, res) => {
            res.json("test router")
        })

    return router;
}

const get_institutions = () : Router => {
    let router = express.Router();
    router
        .get("/institutions", async (req: any, res) => {
            try {
                const institutions = await institutionModel.find({ "memberIds._id" : req.user_id});
                return res.status(200).json(institutions!);
            } catch (err) {
                console.error(`Error getting institutions for user ${req.user_id}`);
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to get institutions");
            }
        })

    return router;
}

const get_institution = () : Router => {
    let router = express.Router();
    router
        .get("/institution/:id", check_auth(), async (req: any, res) =>{
            const institutionId = req.params.id;
            try {
                const institution = await institutionModel.findById(institutionId);
                return res.status(200).json(institution);
            } catch (err) {
                console.error(JSON.stringify(err));
                return res.status(404).send(`Institution ${institutionId} not found`);
            }
        })
    return router;
}

export { create_institution, test_institution, get_institutions, get_institution }
