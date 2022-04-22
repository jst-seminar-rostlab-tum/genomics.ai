import express, {Router} from "express";
import ProjectService from "../../../../database/services/project.service";
import { AddProjectDTO } from "../../../../database/dtos/project.dto";
import UserService from "../../../../database/services/user.service";
import {ObjectId} from "mongoose";
import { visibilityStatus } from "../../../../database/models/project";
import check_auth from "../../middleware/check_auth";
import {mailer} from "../../../../util/mailer";

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
            // // Is it not possible that there exists projects with same names?

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
                return res.status(500).send("Unable to create the project.");
            }
        })

    return router;
}

const invite_person_to_a_project = (): Router => {
    let router = express.Router();

    router.put("/projects/:id/invite", check_auth(), async (req: any, res) => {
        try {
            const {userId}: {userId: ObjectId}  = req.body;
            const projectId: string = req.params.id;

            if(!(userId && projectId))
                return res.status(400).send("Missing parameters.");

            const user = await UserService.getUserById(userId);
            if (! user )
                return res.status(400).send("User to be invited does not exist.");

            const project = await ProjectService.getProjectById(projectId);
            if (! project)
                return res.status(400).send("Project does not exist.");

            var tempUserId = String(userId);
            var tempListAdmins = project.adminIds.map(String);
            var tempListMembers = project.memberIds.map(String);

            if (tempListMembers.includes(tempUserId))
                return res.status(409).send("User is already a member of the project")

            if (tempListAdmins.includes(tempUserId))
                return res.status(409).send("User is an admin of the project")

            try {
                const project_updated = await ProjectService.addInvitationMemberId(projectId, userId);
                if (!project_updated)
                    return res.status(400).send("Error when adding the user to members of the project.");

                try {
                    await mailer.send(user.email, "[GeneCruncher] Invitation to a project", "invitation_to_project", {
                        firstname: user.firstName,
                        projectname: project.title
                        })
                } catch (e) {
                    console.error("Error when sending invitation of user to a project.")
                    console.error(JSON.stringify(e));
                    console.error(e);
                    return res.status(500).send("Error when sending email. Invitation has been stored.");
                }

                return res.status(200).json("Invitation has been sent successfully.");

            } catch(err) {
                console.error("Error when trying to register invitation of user to a project.")
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

const add_user_to_admin = (): Router => {
    let router = express.Router();

    router.put("/projects/:id/admin", check_auth(), async (req: any, res) => {
        try {
            const {userId}: {userId: ObjectId}  = req.body;
            const projectId: string = req.params.id;

            if(!(userId && projectId))
                return res.status(400).send("Missing parameters.");

            const user = await UserService.getUserById(userId)
            if (! user )
                return res.status(400).send("User does not exist.");

            const project = await ProjectService.getProjectById(projectId);
            if (! project)
                return res.status(400).send("Project does not exist.");

            var tempUserId = String(userId);
            var tempListAdmins = project.adminIds.map(String);
            var tempListMembers = project.memberIds.map(String);

            if (tempListAdmins.includes(tempUserId))
                return res.status(409).send("User is already an admin.")

            if (!tempListMembers.includes(tempUserId))
                return res.status(409).send("User is not a member of the project. It should be first a member to become an admin.")            

            try {
                const project_updated = await ProjectService.addAdminToProject(projectId, userId);

                if (!project_updated)
                    return res.status(400).send("Error when changing to user to admin profile.");

                return res.status(200).json("User has been changed to admin.");

            } catch(err) {
                console.error("Error when trying to register new admin of a given a project.")
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
}

const join_member = (): Router => {
    let router = express.Router();

    router.put("/projects/:id/join", check_auth(), async (req: any, res) => {
        try {
            const {userId}: {userId: ObjectId}  = req.body;
            const projectId: string = req.params.id;

            if(!(userId && projectId))
                return res.status(400).send("Missing parameters.");

            const user = await UserService.getUserById(userId)
            if (! user )
                return res.status(400).send("User does not exist.");
            if ( !user.isEmailVerified )
                return res.status(409).send("User has not been verified.")

            const project = await ProjectService.getProjectById(projectId);
            if (! project)
                return res.status(400).send("Project does not exist.");

            var tempUserId = String(userId);
            var tempListAdmins = project.adminIds.map(String);
            var tempListMembers = project.memberIds.map(String);

            if (tempListAdmins.includes(tempUserId))
                return res.status(409).send("User is an admin of the project.")

            if (tempListMembers.includes(tempUserId))
                return res.status(409).send("User is already a member of the project.")   
                
            //[PENDING] Consider whether is necessary to validate if a user should have been invited before it is joined. i.e. it should exist a record in invitedMemberIds

            try {
                const project_updated = await ProjectService.addNewMemberIntoProject(projectId, userId);

                if (!project_updated)
                    return res.status(400).send("Error when joining a new member into the project.");

                return res.status(200).json("User has been joined.");

            } catch(err) {
                console.error("Error when trying to join a new member into a given project.")
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
}

export { create_project, invite_person_to_a_project, add_user_to_admin, join_member }
