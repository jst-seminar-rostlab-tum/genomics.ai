import express, {Router} from "express";
import {ExtRequest} from "../../../../definitions/ext_request";
import ProjectService from "../../../../database/services/project.service";
import TeamService from "../../../../database/services/team.service";
import InstitutionService from "../../../../database/services/institution.service";
import check_auth from "../../middleware/check_auth";
import {ObjectId} from "mongoose";

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

/**
 *  Get a specific project by its id
 */
const get_project_by_id = (): Router => {
    let router = express.Router();

    router.get("/project/:id", check_auth(), async (req: ExtRequest, res: any) => {
        try {
            const projectId: string = req.params.id;

            if(!projectId)
                return res.status(400).send("Missing parameters.");

            const project = await ProjectService.getProjectById(projectId);
            if (!project)
                return res.status(400).send("Project does not exist.");

            let allowAccess = false;

            /* is user the owner of the project */
            let isOwner: boolean = req.user_id! === await ProjectService.getOwner(projectId);
            allowAccess ||= isOwner;

            /* is the project part of a team */
            let teamId: (ObjectId | null) = await ProjectService.getTeam(projectId);

            if (!allowAccess && teamId) {
              const visibility: (string | null) = await TeamService.getVisibility(teamId);
              const isAdmin = await TeamService.isAdmin(req.user_id!, teamId);
              const isMember = await TeamService.isMember(req.user_id!, teamId);

              /* is the project part of a PUBLIC team */
              const inPublicTeam: boolean = "PUBLIC" === visibility;

              /* is the user part of the team that the project belongs to */
              const userInTeam: boolean = isAdmin || isMember;

              if ("BY_INSTITUTION" === visibility) {
                  /* is the project part of a team that is a PUBLIC institution */
                  const institutionId = await TeamService.getInstitution(teamId);
                  const institutionVisibility = await InstitutionService.getVisibility(institutionId!);
                  const inPublicInstitution: boolean = "PUBLIC" === institutionVisibility;

                  /* is the user part of the institution */
                  const isAdmin = await InstitutionService.isAdmin(req.user_id!, institutionId!);
                  const isMember = await InstitutionService.isMember(req.user_id!, institutionId!);
                  const userInInstitution: boolean = isAdmin || isMember;

                  allowAccess ||= inPublicInstitution || userInInstitution;
              }

              allowAccess ||= inPublicTeam || userInTeam;
            }

            if (!allowAccess)
                return res.status(401).send("Access Denied");

            return res.status(200).json(project);

        } catch(e) {
            console.error("Error in get_project_by_id()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error. Cannot get project by id.");
        }
    })
    return router;
};

export { get_userProjects, get_project_by_id };
