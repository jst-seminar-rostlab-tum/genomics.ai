import express, {Router} from "express";
import ProjectService from "../../../../database/services/project.service";
import TeamService from "../../../../database/services/team.service";
import { AddTeamDTO } from "../../../../database/dtos/team.dto";
import UserService from "../../../../database/services/user.service";
import InstitutionService from "../../../../database/services/institution.service";
import {ObjectId} from "mongoose";
import { visibilityStatus } from "../../../../database/models/team";
import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";
import {mailer} from "../../../../util/mailer";

const create_team = () : Router => {
    let router = express.Router();

    router
        .post("/teams", check_auth(), async (req: any, res) => {

            const {title, description, visibility} = req.body;
            const admin_user_id = req.user_id;

            if (!(title && description && visibility && admin_user_id))
                return res.status(400).send("Missing parameters");

            if(!(Object.values(visibilityStatus).includes(visibility)) )
                return res.status(400).send("Visibility parameter is wrong format. You should type one of the followings: PRIVATE, PUBLIC, BY_INSTITUTION");

            // TODO: make this check later
            // const team = await TeamService.getTeamByTitle(title);
            // if (team)
            //     return res.status(409).send("Team with the given name already exists!");
            // // Is it not possible that there exists teams with same names?

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
};

/**
 *  Invites a user to a team by adding this userId to
 *  the invitedMemberIds of this team.
 */
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

            const isAdmin: boolean = await TeamService.isAdmin(userId, team);
            const isMember: boolean = await TeamService.isMember(userId, team);

            if (isAdmin || isMember)
                return res.status(409).send("User is already part of the team")

            try {
                const team_updated = await TeamService.addInvitationMemberId(teamId, userId);
                if (!team_updated)
                    return res.status(400).send("Error when adding the user to members of the team.");

                try {
                    await mailer.send(user.email, "[GeneCruncher] Invitation to a team", "invitation_to_team", {
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
};

/**
 *  Adds the given projectId to the projects list of the team.
 */
const add_project_to_team = (): Router => {
    let router = express.Router();

    router.put("/teams/:id/add_project", check_auth(), async (req: ExtRequest, res: any) => {
        try {
            const {projectId}: {projectId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(projectId && teamId))
                return res.status(400).send("Missing parameters.");

            const project = await ProjectService.getProjectById(projectId);
            if (!project)
                return res.status(400).send("Project to be added does not exist.");

            const team = await TeamService.getTeamById(teamId);
            if (!team)
                return res.status(400).send("Team does not exist.");

            /* the user should be either an admin or a member of the team */
            const isAdmin: boolean = await TeamService.isAdmin(req.user_id!, team);
            const isMember: boolean = await TeamService.isMember(req.user_id!, team);
            if (!(isAdmin || isMember))
                return res.status(401).send("Unauthenticated User! The user is not part of the team.");

            if (project.owner!.toString() !== req.user_id!.toString()) // cannot be null since the project must exist at this point
                return res.status(401).send("User is not the project owner!");

            try {
                const team_updated = await TeamService.addProject(teamId, projectId);
                const project_updated = await ProjectService.setTeamOfProject(projectId, teamId);

                if (!team_updated)
                    return res.status(500).send("Error when adding the project to project list of the team.");

                if (!project_updated)
                    return res.status(500).send("Error when setting the team id of the project.");

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
};

/**
 *  Add the given userId to the admin list of a team.
 */
const add_user_to_admin = (): Router => {
    let router = express.Router();

    router.put("/teams/:id/admin", check_auth(), async (req: any, res) => {
        try {
            const {userId}: {userId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(userId && teamId))
                return res.status(400).send("Missing parameters.");

            const user = await UserService.getUserById(userId)
            if (! user )
                return res.status(400).send("User does not exist.");

            const team = await TeamService.getTeamById(teamId);
            if (! team)
                return res.status(400).send("Team does not exist.");

            const isAdmin: boolean = await TeamService.isAdmin(userId, team);
            const isMember: boolean = await TeamService.isMember(userId, team);
            if (isAdmin)
                return res.status(409).send("User is already an admin.")
            if (!isMember)
                return res.status(409).send("User is not a member of the team. It should first become a member to become an admin.")

            try {
                const team_updated = await TeamService.addAdminToTeam(teamId, userId);

                if (!team_updated)
                    return res.status(400).send("Error when changing to user to admin profile.");

                return res.status(200).json("User has been changed to admin.");

            } catch(err) {
                console.error("Error when trying to register new admin of a given a team.")
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to register new admin.");
            }
        } catch(e) {
            /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
            console.error("Error in add_user_to_admin()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }
    })

    return router;
};

/**
 *  Add the given userId to the member list of a team.
 */
const join_member = (): Router => {
    let router = express.Router();

    router.put("/teams/:id/join", check_auth(), async (req: any, res) => {
        try {
            const {userId}: {userId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(userId && teamId))
                return res.status(400).send("Missing parameters.");

            const user = await UserService.getUserById(userId)
            if (! user )
                return res.status(400).send("User does not exist.");
            if ( !user.isEmailVerified )
                return res.status(409).send("User has not been verified.")

            const team = await TeamService.getTeamById(teamId);
            if (! team)
                return res.status(400).send("Team does not exist.");

            const isAdmin = await TeamService.isAdmin(userId, team);
            if (isAdmin)
                return res.status(409).send("User is an admin of the team.")

            const isMember = await TeamService.isMember(userId, team);
            if (isMember)
                return res.status(409).send("User is already a member of the team.")

            //[PENDING] Consider whether is necessary to validate if a user should have been invited before it is joined. i.e. it should exist a record in invitedMemberIds

            try {
                const team_updated = await TeamService.addNewMemberIntoTeam(teamId, userId);

                if (!team_updated)
                    return res.status(400).send("Error when joining a new member into the team.");

                return res.status(200).json("User has been joined.");

            } catch(err) {
                console.error("Error when trying to join a new member into a given team.")
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to register new admin.");
            }
        } catch(e) {
            /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
            console.error("Error in join_member()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }

    })

    return router;
};

/**
 *  Add a team to an institution
 */
const add_team_to_institution = (): Router => {
    let router = express.Router();

    router.put("/teams/:id/institution", check_auth(), async (req: any, res) => {
        try {
            const {institutionId}: {institutionId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(institutionId && teamId))
                return res.status(400).send("Missing parameters.");

            const institution = await InstitutionService.getInstitutionById(institutionId)
            if (! institution )
                return res.status(400).send("Institution does not exist.");

            const team = await TeamService.getTeamById(teamId);
            if (!team)
                return res.status(400).send("Team does not exist.");

            if (team.institutionId)
                return res.status(409).send("Team has been already associated with an institution.")

            try {
                const institution_updated = await TeamService.setInstitutionOfTeam(teamId, institutionId);

                if (!institution_updated)
                    return res.status(400).send("Error when associating the team with the institution.");

                const team2 = await TeamService.getTeamById(teamId);
                return res.status(200).json(team2);

            } catch(err) {
                console.error("Error when trying to join the team into the institution.")
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to associate the team with the institution.");
            }
        } catch(e) {
            /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
            console.error("Error in add_team_to_institution()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }

    })

    return router;
};

/**
 *  Remove a team from an institution
 */
const remove_team_from_institution = (): Router => {
    let router = express.Router();

    router.delete("/teams/:id/institution", check_auth(), async (req: any, res) => {
        try {
            const {institutionId}: {institutionId: ObjectId}  = req.body;
            const teamId: string = req.params.id;

            if(!(institutionId && teamId))
                return res.status(400).send("Missing parameters.");

            const institution = await InstitutionService.getInstitutionById(institutionId)
            if (! institution )
                return res.status(400).send("Institution does not exist.");

            const team = await TeamService.getTeamById(teamId);
            if (! team)
                return res.status(400).send("Team does not exist.");

            try {
                const institution_updated = await InstitutionService.removeTeamFromInstitution(teamId, institutionId);

                if (!institution_updated)
                    return res.status(400).send("Error when removing the team from the institution.");

                return res.status(200).json("Team has been removed.");

            } catch(err) {
                console.error("Error when trying to join the team into the institution.")
                console.error(JSON.stringify(err));
                console.error(err);
                return res.status(500).send("Unable to register the team into the institution.");
            }
        } catch(e) {
            /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
            console.error("Error in add_team_to_institution()")
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Internal error.");
        }

    })

    return router;
};

const get_teams = () : Router => {
    let router = express.Router();
    router
        .get("/teams", check_auth(), async (req: any, res) => {
            const query = { ...req.query };
            try{
                const teams = await TeamService.getTeams(query);
                return res.status(200).json(teams);
            } catch (err){
                console.error(JSON.stringify(err));
                return res.status(404).send(`No teams found`);
            }
        })

    return router;
}

export { create_team, invite_person_to_a_team, add_user_to_admin, join_member, add_team_to_institution, remove_team_from_institution, add_project_to_team, get_teams};
