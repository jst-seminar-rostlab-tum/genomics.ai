import express, {Router} from "express";
import check_auth from "../../middleware/check_auth";
import ProjectService from "../../../../database/services/project.service";

const get_projects = () : Router => {
    let router = express.Router();
    router
        .get("/projects", check_auth(), async( req: any, res) => {
            const query = { ...req.query };
            try{
                const projects = await ProjectService.getProjects(query);
                return res.status(200).json(projects);
            } catch (err){
                console.error(JSON.stringify(err));
                return res.status(404).send(`No project found`);
            }
        })
    return router;
}

export { get_projects }