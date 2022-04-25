import express, {Router} from "express";
import {ExtRequest} from "../../../../definitions/ext_request";
import ProjectService from "../../../../database/services/project.service";
import check_auth from "../../middleware/check_auth";

/**
 *  Returns all the projects of the user in order of upload date.
 *  @param sort? - query parameter, specificies the order of the sort
 */
const get_userProjects = () : Router => {
    let router = express.Router();

    router.get("/projects", check_auth(), async (req: ExtRequest, res: any) => {
        try {
            const { sort, ...rest } = req.params;
            let sortParam: number;
            if (sort === '1')
                sortParam = 1;
            else
                sortParam = -1;
            const projects = req.user_id === undefined ? [] :
              await ProjectService.getProjectByOwner(req.user_id, sortParam);
            return res.status(200).json(projects);
        } catch (e) {
            console.error("Error in get_userProjects_route()");
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }
    });
    return router;
}

export { get_userProjects };
