import express, {Router} from "express";
import ProjectService from "../../../../database/services/project.service";
import { AddProjectDTO } from "../../../../database/dtos/project.dto";
import UserService from "../../../../database/services/user.service";
import { visibilityStatus } from "../../../../database/models/project";
import check_auth from "../../middleware/check_auth";

const create_project = () : Router => {
    let router = express.Router();

    router
        .post("/projects", check_auth(), async (req: any, res) => {

            const {title, description, visibility} = req.body;
            const admin_user_id = req.user_id;

            if (!(title && description && visibility && admin_user_id))
                return res.status(400).send("Missing parameters");

            if(!(Object.values(visibilityStatus).includes(visibility)) )
                return res.status(400).send("Visibility parameter is wrong format. You should type one of the followings: PRIVATE, PUBLIC, BY_INSTITUTION");

            // TODO: make this check later
            // const project = await ProjectService.getProjectByTitle(title);
            // if (project)
            //     return res.status(409).send("Project with the given name already exists!");

            const admin = await UserService.getUserById(admin_user_id);
            if (!admin)
                return res.status(404).send("Admin that you are trying to assign does not exists!");

            try{
                const projectToAdd: AddProjectDTO = {
                    title,
                    description,
                    visibility,
                    adminIds: [admin_user_id]
                };
                const project = await ProjectService.addProject(projectToAdd);

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

export { create_project }
