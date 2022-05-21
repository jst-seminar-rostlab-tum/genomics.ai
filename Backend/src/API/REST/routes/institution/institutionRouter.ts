import express, { Router, Request, Response } from "express";

import { Schema, ObjectId } from "mongoose";

import InstitutionService from "../../../../database/services/institution.service";
import TeamsService from "../../../../database/services/team.service";
import { AddInstitutionDTO, UpdateInstitutionDTO } from "../../../../database/dtos/institution.dto";
import UserService from "../../../../database/services/user.service";

import check_auth from "../../middleware/check_auth";
import { validationMdw } from "../../middleware/validation";
import { institution_admin_auth } from "../../middleware/check_institution_auth";
import { mailer } from "../../../../util/mailer";

const create_institution = (): Router => {
  let router = express.Router();

  router.post("/institutions", validationMdw, check_auth(), async (req: any, res) => {
    const { name, country, profilePictureURL, backgroundPictureURL } = req.body;
    const admin_user_id = req.user_id;

    if (!(name && country && admin_user_id)) return res.status(400).send("Missing parameters");

    const institution = await InstitutionService.getInstitutionByName(name);
    if (institution) return res.status(409).send("Institution with the given name already exists!");

    const user = await UserService.getUserById(admin_user_id);
    if (!user) return res.status(404).send("Admin that you are trying to assign does not exists!");

    try {
      const institutionToAdd: AddInstitutionDTO = {
        name,
        country,
        profilePictureURL,
        backgroundPictureURL,
        adminIds: [admin_user_id],
      };
      const institution = await InstitutionService.addInstitution(institutionToAdd);

      return res.status(201).json(institution);
    } catch (err) {
      console.error("Error registering institution!");
      console.error(JSON.stringify(err));
      console.error(err);
      return res.status(500).send("Unable to create institution. (DB-error)");
    }
  });

  return router;
};

const update_institution = (): Router => {
  let router = express.Router();

  router.put(
    "/institutions/:id",
    validationMdw,
    check_auth(),
    institution_admin_auth,
    async (req: any, res) => {
      const { name, country } = req.body;
      const institution_to_be_updated_id = req.params.id;

      if (!(name || country)) return res.status(400).send("Missing parameters");

      try {
        const institution = await InstitutionService.getInstitutionByName(name);
        if (institution)
          return res.status(409).send("Institution with the given name already exists!");

        const institutionToUpdate: UpdateInstitutionDTO = {
          name,
          country,
        };
        await InstitutionService.updateInstitution(
          institution_to_be_updated_id,
          institutionToUpdate
        );

        const updatedInstitution = await InstitutionService.getInstitutionById(
          institution_to_be_updated_id
        );
        return res.status(200).send(updatedInstitution); 
      } catch (err) {
        console.error("Error updating institution!");
        console.error(JSON.stringify(err));
        console.error(err);
        return res.status(500).send("Unable to update institution. (DB-error)");
      }
    }
  );

  return router;
};

const invite_to_institution = (): Router => {
  let router = express.Router();

  router.put(
    "/institutions/:id/invite",
    validationMdw,
    check_auth(),
    async (req: any, res: any) => {
      const { userId, email }: { userId: Schema.Types.ObjectId; email: string } = req.body;
      const institutionId_to_modify = req.params.id;

      try {
        if (!(institutionId_to_modify && (userId || email)))
          return res.status(400).send("Missing parameters.");

        var user;
        if (userId) {
          user = await UserService.getUserById(userId);
          console.log(user);
        } else {
          user = await UserService.getUserByEmail(email, false);
          console.log(user);
        }
        console.log(user);

        if (!user) return res.status(404).send("User to be invited does not exist.");

        if (await InstitutionService.findMemeberOrInvitedById(userId, institutionId_to_modify))
          return res
            .status(404)
            .send(
              "User that you are trying to invite to this institution already is an invited member or is a member!"
            );

        const updatedInstitution = await InstitutionService.inviteToInstitution(
          institutionId_to_modify,
          userId
        );

        if (updatedInstitution) {
          try {
            await mailer.send(
              user.email,
              "[GeneCruncher] Invitation to an institution",
              "invitation_to_institution",
              {
                institution: updatedInstitution.name,
                country: updatedInstitution.country,
                firstname: user.firstName,
              }
            );
          } catch (e) {
            console.error("Error when sending invitation of user to an institution.");
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Error when sending email. Invitation has been stored.");
          }
          return res.json(updatedInstitution);
        } else {
          return res.status(409).send("Could not invite person to institution!");
        }
      } catch (e) {
        return res.status(500).send("Error: Something went wrong internal error!");
      }
    }
  );

  return router;
};

const make_user_admin_of_institution = (): Router => {
  let router = express.Router();

  router.put(
    "/institutions/:id/admin",
    validationMdw,
    check_auth(),
    institution_admin_auth,
    async (req: any, res) => {
      const { userId }: { userId: Schema.Types.ObjectId } = req.body;
      const institutionId_to_modify = req.params.id;

      try {
        if (!userId) return res.status(400).send("Missing parameter");

        if (!(await UserService.getUserById(userId)))
          return res.status(404).send("User that you are trying to make as admin does not exists!");

        const institutionToBeUpdated = await InstitutionService.findMemeberById(
          userId,
          institutionId_to_modify
        );

        if (!institutionToBeUpdated)
          return res.status(409).send("User that you are trying to make as admin is not a member!");

        if (institutionToBeUpdated?.adminIds.includes(userId))
          return res.status(409).send("User is already an admin!");

        const updatedInstitution = await InstitutionService.makeUserAnAdminOfInstitution(
          institutionId_to_modify,
          userId
        );

        if (updatedInstitution) {
          res.json(updatedInstitution);
        } else {
          return res.status(409).send("Could not make the user admin of the institution!");
        }
      } catch (error) {
        return res.status(500).send("Something went wrong: " + error);
      }
    }
  );

  return router;
};

const join_as_member_of_institution = (): Router => {
  let router = express.Router();

  router.put("/institutions/:id/join", check_auth(), async (req: any, res) => {
    const institutionId_to_modify = req.params.id;
    const current_user = req.user_id;
    try {
      const institutionToBeUpdated = await InstitutionService.getInstitutionById(
        institutionId_to_modify
      );

      if (!institutionToBeUpdated?.invitedMemberIds.includes(current_user))
        return res.status(409).send("Could not join as you are not an invited member!");

      const updatedInstitution = await InstitutionService.makeUserMemberOfInstitution(
        institutionId_to_modify,
        current_user
      );

      if (updatedInstitution) {
        res.json(updatedInstitution);
      } else {
        return res.status(409).send("Could not join as member of the institution!");
      }
    } catch (error) {
      return res.status(500).send("Something went wrong: " + error);
    }
  });

  return router;
};

const get_institution = (): Router => {
  let router = express.Router();
  router.get("/institutions/:id", check_auth(), async (req: any, res) => {
    const institutionId = req.params.id;
    try {
      const institution = await InstitutionService.getInstitutionById(institutionId);

      if (institution != null) return res.status(200).json(institution);
      return res.status(404).send(`Institution ${institutionId} not found`);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).send(`Internal server error`);
    }
  });
  return router;
};

const get_institutions = (): Router => {
  let router = express.Router();
  router.get("/institutions", check_auth(), async (req: any, res) => {
    const query = { ...req.query };
    try {
      const institutions = await InstitutionService.filterInstitutions(query);

      if (institutions != null) return res.status(200).json(institutions);
      return res.status(404).send(`No institutions found`);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).send(`Internal server error`);
    }
  });
  return router;
};

const get_members_of_institution = (): Router => {
  let router = express.Router();
  router.get("/institutions/:id/members", check_auth(), async (req: Request, res: Response) => {
    const institutionId = req.params.id;
    try {
      const institution = await InstitutionService.getMembersOfInstitution(institutionId);
      if (institution == null) {
        return res.status(404).send(`Institution ${institutionId} not found`);
      }
      return res.status(200).send(institution.memberIds);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).json({ error: "General server error" });
    }
  });
  return router;
};

const get_teams_of_institution = (): Router => {
  let router = express.Router();
  router.get("/institutions/:id/teams", check_auth(), async (req: Request, res: Response) => {
    const institutionId = req.params.id;
    try {
      const teams = await TeamsService.getTeams({ institutionId: institutionId });

      return res.status(200).send(teams);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).json({ error: "General server error" });
    }
  });
  return router;
};

const get_projects_of_institution = (): Router => {
  let router = express.Router();
  router.get("/institutions/:id/projects", check_auth(), async (req: Request, res: Response) => {
    const institutionId = req.params.id;
    try {
      const projects = await InstitutionService.getProjectsOfInstitution(institutionId);
      if (projects == null) {
        return res.status(404).send(`Institution ${institutionId} not found`);
      }
      return res.status(200).send(projects);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).json({ error: "General server error" });
    }
  });
  return router;
};
const get_users_institutions = (): Router => {
  let router = express.Router();
  router.get("/users/:id/institutions", check_auth(), async (req: any, res) => {
    const userId = req.params.id;
    try {
      const institutions = await InstitutionService.getUsersInstitutions(userId);
      if (institutions != null) return res.status(200).json(institutions);
      return res.status(404).send(`No institutions found`);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).send(`Internal server error`);
    }
  });
  return router;
};

const disjoin_member_of_institution = (): Router => {
  let router = express.Router();

  router.delete("/institutions/:id/join", validationMdw, check_auth(), async (req: any, res) => {
    try {
      const { userId }: { userId: ObjectId } = req.body;
      const institutionId: string = req.params.id;
      const user_id_jwt = req.user_id;

      if (!(userId && institutionId)) return res.status(400).send("Missing parameters.");

      const user = await UserService.getUserById(userId);
      if (!user) return res.status(409).send("User does not exist.");
      if (userId != user_id_jwt)
        return res.status(409).send("Information of the user does not match.");

      const institution = await InstitutionService.getInstitutionById(institutionId);
      if (!institution) return res.status(409).send("Institution does not exist.");

      var tempUserId = String(userId);
      var tempListAdmins = institution.adminIds.map(String);
      var tempListMembers = institution.memberIds.map(String);

      if (!(tempListMembers.includes(tempUserId) || tempListAdmins.includes(tempUserId)))
        return res.status(409).send("You are not member of the institution.");

      if (tempListAdmins.includes(tempUserId) && tempListAdmins.length == 1)
        return res.status(403).send("You are the only one admin of the institution.");

      try {
        const inst_updated = await InstitutionService.removeMemberFromTeam(institutionId, userId);

        if (!inst_updated)
          return res.status(500).send("Error when removing you from the institution.");

        const instRes = await InstitutionService.getInstitutionById(institutionId);
        instRes.memberIds.push(...instRes.adminIds);

        return res.status(200).json(instRes);
      } catch (err) {
        console.error("Error when trying to remove a member from a institution.");
        console.error(JSON.stringify(err));
        console.error(err);
        return res.status(500).send("Unable to remove user from a institution.");
      }
    } catch (e) {
      /* Added since a test proved that if user sends a request with incorrect parameter names, it is able to shutdown the server. */
      console.error("Error in disjoin_member_of_institution()");
      console.error(JSON.stringify(e));
      console.error(e);
      return res.status(500).send("Internal error.");
    }
  });

  return router;
};

export {
  create_institution,
  update_institution,
  invite_to_institution,
  make_user_admin_of_institution,
  join_as_member_of_institution,
  get_institutions,
  get_institution,
  get_members_of_institution,
  get_teams_of_institution,
  get_projects_of_institution,
  get_users_institutions,
  disjoin_member_of_institution,
};
