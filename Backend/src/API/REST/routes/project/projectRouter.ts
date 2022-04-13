import express, {Router} from "express";
import { IProject, projectModel } from "../../../../database/models/project";
import { userModel } from "../../../../database/models/user";
import { visibilityStatus } from "../../../../database/models/project";

export default function project_route() : Router {
    let router = express.Router();

    router
        .post("/projects", async (req: any, res) => {
            
            const {title, description, adminId, invitedMemberIds, memberIds, visibility, projectJobs} = req.body;

            if (!(title && description && visibility && adminId))
                return res.status(400).send("Missing parameters");

            if (await projectModel.findOne({title}))
                return res.status(409).send("Project with the given name already exists!");

            if (! await userModel.findOne({_id: adminId}))
                return res.status(404).send("Admin that you are trying to assign does not exists!");   

            let project : (IProject | undefined) = undefined;

            // TODO: think again the way of getting ids response and how to store the ids 
            try{
                project = await projectModel.create({
                    title,
                    description,
                    invitedMemberIds,
                    memberIds,
                    visibility,
                    projectJobs,
                    adminIds: [adminId]
                });

                return res.status(201).json(project);
            }catch(err){
                console.error("Error registering project!");
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to create project. (DB-error)");
            }
        })

    return router;
}
