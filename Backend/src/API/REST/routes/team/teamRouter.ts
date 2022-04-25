import express, {Router} from "express";
import ProjectService from "../../../../database/services/project.service";
import TeamService from "../../../../database/services/team.service";
import { AddTeamDTO } from "../../../../database/dtos/team.dto";
import UserService from "../../../../database/services/user.service";
import {ObjectId} from "mongoose";
import { visibilityStatus } from "../../../../database/models/team";
import check_auth from "../../middleware/check_auth";
import {mailer} from "../../../../util/mailer";

const create_team = () : Router => {
    let router = express.Router();

    router
        .post("/teams", async (req: any, res) => {

            const {title, description, visibility} = req.body;
            const admin_user_id = req.user_id;

            if (!(title && description && visibility && admin_user_id))
                return res.status(400).send("Missing parameters");

            if(!(Object.values(visibilityStatus).includes(visibility)) )
                return res.status(400).send("Visibility parameter is wrong format. You should type one of the followings: PRIVATE, PUBLIC, BY_INSTITUTION");

            // TODO: make this check later
            // const team = await TeamService.getTeamByTitle(title);
            // if (team)
            //     return res.status(409).send("Project with the given name already exists!");
            // // Is it not possible that there exists projects with same names?

            const admin = await UserService.getUserById(admin_user_id);
            if (!admin)
                return res.status(404).send("Admin that you are trying to assign does not exists!");

            try{
                const teamToAdd: AddTeamDTO = {
                    title,
                    description,
                    visibility,
                    adminIds: [admin_user_id]
                };
                const project = await TeamService.addTeam(teamToAdd);

                return res.status(201).json(project);
            }catch(err){
                console.error("Error registering team!");
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to create the team.");
            }
        })

    return router;
}

const invite_person_to_a_team = (): Router => {
    let router = express.Router();

    router.put("/teams/:id/invite", check_auth(), async (req: any, res) => {
        try {
            const {userId}: {userId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(userId && teamId))
                return res.status(400).send("Missing parameters.");

            const user = await UserService.getUserById(userId);
            if (! user )
                return res.status(400).send("User to be invited does not exist.");

            const team = await TeamService.getTeamById(teamId);
            if (! team)
                return res.status(400).send("Team does not exist.");

            var tempUserId = String(userId);
            var tempListAdmins = team.adminIds.map(String);
            var tempListMembers = team.memberIds.map(String);

            if (tempListMembers.includes(tempUserId))
                return res.status(409).send("User is already a member of the team")

            if (tempListAdmins.includes(tempUserId))
                return res.status(409).send("User is an admin of the team")

            try {
                const team_updated = await TeamService.addInvitationMemberId(teamId, userId);
                if (!team_updated)
                    return res.status(400).send("Error when adding the user to members of the team.");

                try {
                    await mailer.send(user.email, "[GeneCruncher] Invitation to a team", "invitation_to_project", {
                        firstname: user.firstName,
                        teamname: team.title
                        })
                } catch (e) {
                    console.error("Error when sending invitation of user to a team.")
                    console.error(JSON.stringify(e));
                    console.error(e);
                    return res.status(500).send("Error when sending email. Invitation has been stored.");
                }

                return res.status(200).json("Invitation has been sent successfully.");

            } catch(err) {
                console.error("Error when trying to register invitation of user to a team.")
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to send invitation to the desired user.");
            }
        } catch(e) {
            /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
            console.error("Error in invite_person_to_a_project()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }

    })

    return router;
}

const add_project_to_team = (): Router => {
    let router = express.Router();

    router.put("/teams/:id/add_project", check_auth(), async (req: any, res) => {
        try {
            const {projectId}: {projectId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(projectId && teamId))
                return res.status(400).send("Missing parameters.");

            const project = await ProjectService.getProjectById(projectId);
            if (! project )
                return res.status(400).send("Project to be added does not exist.");

            const team = await TeamService.getTeamById(teamId);
            if (! team)
                return res.status(400).send("Team does not exist.");

            try {
                const team_updated = await TeamService.addProject(teamId, projectId);

                if (!team_updated)
                    return res.status(400).send("Error when adding the project to project list of the team.");

                return res.status(200).json("The project is successfully added.");
            } catch(err) {
                console.error("Error when trying to add the project to a team.")
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to add the project to the team.");
            }
        } catch(e) {
            /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
            console.error("Error in invite_person_to_a_project()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }
    })
    return router;
}

export { create_team, invite_person_to_a_team, add_project_to_team };
