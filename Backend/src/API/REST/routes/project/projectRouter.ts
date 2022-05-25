import express, { Router } from "express";
import check_auth from "../../middleware/check_auth";
import { ExtRequest } from "../../../../definitions/ext_request";
import ProjectService from "../../../../database/services/project.service";
import TeamService from "../../../../database/services/team.service";
import InstitutionService from "../../../../database/services/institution.service";
import { ObjectId } from "mongoose";
import { validationMdw } from "../../middleware/validation";
import ProjectUpdateTokenService from "../../../../database/services/project_update_token.service";
import DeletedProjectService from "../../../../database/services/deletedProject.service";
import { IDeletedProject } from "../../../../database/models/deleted_projects";
import { UpdateProjectDTO } from "../../../../database/dtos/project.dto";
import s3 from "../../../../util/s3";
import { DeleteObjectRequest } from "aws-sdk/clients/s3";
import { ProjectStatus } from "../../../../database/models/project";
import { query_path, result_model_path, result_path } from "../file_upload/bucket_filepaths";

const get_projects = (): Router => {
  let router = express.Router();
  router.get("/projects", check_auth(), async (req: any, res) => {
    const query = { ...req.query };
    try {
      const projects = await ProjectService.getProjects(query);

      if (projects != null) return res.status(200).json(projects);
      return res.status(404).send(`No project found`);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).send(`Internal server error`);
    }
  });
  return router;
};

/**
 *  Returns all the projects of the user in order of upload date.
 *  @param sort? - query parameter, specificies the order of the sort
 */
const get_userProjects = (): Router => {
  let router = express.Router();

  router.get("/ownprojects", check_auth(), async (req: ExtRequest, res: any) => {
    try {
      const { sort, ...rest } = req.params;
      let sortParam: number;
      if (sort === "1") sortParam = 1;
      else sortParam = -1;
      const projects =
        req.user_id === undefined
          ? []
          : await ProjectService.getProjectByOwner(req.user_id, sortParam);
      return res.status(200).json(projects);
    } catch (e) {
      console.error("Error in get_userProjects_route()");
      console.error(JSON.stringify(e));
      console.error(e);
      return res.status(500).send("Internal error.");
    }
  });
  return router;
};

/**
 *  Get a specific project by its id
 */
const get_project_by_id = (): Router => {
  let router = express.Router();

  router.get("/project/:id", check_auth(), async (req: ExtRequest, res: any) => {
    try {
      const projectId: string = req.params.id;

      if (!projectId) return res.status(400).send("Missing parameters.");

      const project = await ProjectService.getProjectById(projectId);
      if (!project) return res.status(400).send("Project does not exist.");

      let allowAccess = false;

      /* is user the owner of the project */
      let isOwner: boolean = req.user_id!.toString() === project.owner.toString();
      allowAccess ||= isOwner;

      /* is the project part of a team */
      let teamId: ObjectId | undefined = project.teamId;
      const team = await TeamService.getTeamById(teamId);

      if (!allowAccess && team) {
        const visibility: string = team.visibility; // "PRIVATE", "PUBLIC", "BY_INSTITUTION"
        const isAdmin = await TeamService.isAdmin(req.user_id!, team);
        const isMember = await TeamService.isMember(req.user_id!, team);

        /* is the project part of a PUBLIC team */
        const inPublicTeam: boolean = "PUBLIC" === visibility;

        /* is the user part of the team that the project belongs to */
        const userInTeam: boolean = isAdmin || isMember;

        if ("BY_INSTITUTION" === visibility) {
          /* is the project part of a team that is a PUBLIC institution */
          const institution = await InstitutionService.getInstitutionById(team.institutionId!);
          /* since visibility equals BY_INSTITUTION this should always hold actually */
          if (institution) {
            const institutionVisibility = institution.visibility!;
            const inPublicInstitution: boolean = "PUBLIC" === institutionVisibility;

            /* is the user part of the institution */
            const isAdmin = await InstitutionService.isAdmin(req.user_id!, institution!);
            const isMember = await InstitutionService.isMember(req.user_id!, institution!);
            const userInInstitution: boolean = isAdmin || isMember;

            allowAccess ||= inPublicInstitution || userInInstitution;
          }
        }

        allowAccess ||= inPublicTeam || userInTeam;
      }

      if (!allowAccess) return res.status(401).send("Access Denied");

      return res.status(200).json(project);
    } catch (e) {
      console.error("Error in get_project_by_id()");
      console.error(JSON.stringify(e));
      console.error(e);
      return res.status(500).send("Internal error. Cannot get project by id.");
    }
  });
  return router;
};

const get_users_projects = (): Router => {
  let router = express.Router();
  router.get("/users/:id/projects", check_auth(), async (req: any, res) => {
    const userId = req.params.id;
    try {
      const projects = await ProjectService.getProjectByOwner(userId);

      if (projects != null) return res.status(200).json(projects);
      return res.status(404).send(`No projects found`);
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).send(`Internal server error`);
    }
  });
  return router;
};

const update_project_results = (): Router => {
  let router = express.Router();
  router.post("/projects/updateresults/:token", validationMdw, async (req, res) => {
    try {
      const updateToken = req.params.token;
      let tokenObject = await ProjectUpdateTokenService.getTokenByToken(updateToken);
      let project = await ProjectService.getProjectById(tokenObject._projectId);
      if (project.status === ProjectStatus.PROCESSING_PENDING) {
        let params: any = {
          Bucket: process.env.S3_BUCKET_NAME!,
          Key: result_path(project._id),
          Expires: 60 * 60 * 24 * 7 - 1, // one week minus one second
        };
        let presignedUrl = await s3.getSignedUrlPromise("getObject", params);
        const updateLocationAndStatus: UpdateProjectDTO = {
          location: presignedUrl,
          status: ProjectStatus.DONE,
        };
        await ProjectService.updateProjectById(project._id, updateLocationAndStatus);
      } else {
        console.log(`Trying to update project with token, but status is already ${project.status}`)
      }
      return res.status(200).send("OK");
    } catch (e) {
      console.error(e);
      return res.status(500).send("Internal server error");
    }
  });
  return router;
};

const delete_project = (): Router => {
  let router = express.Router();
  router.delete("/project/:id", check_auth(), validationMdw, async (req, res) => {
    try {
      const projectId = req.params.id;
      const project = (await ProjectService.getProjectById(projectId))?.toObject();
      if (project == null) return res.status(404).send("Project not found");

      const deletedProject = {
        ...project,
        deletedAt: new Date(),
      };

      await DeletedProjectService.addDeletedProject(deletedProject);
      await ProjectService.deleteProjectById(projectId);

      return res.status(200).send("OK");
    } catch (err) {
      console.error(JSON.stringify(err));
      return res.status(500).send("Internal server error");
    }
  });
  return router;
};

const get_deleted_projects = (): Router => {
  let router = express.Router();
  router.get("/deletedprojects", check_auth(), validationMdw, async (req: ExtRequest, res) => {
    try {
      let projects = await DeletedProjectService.getDeletedProjectsByOwner(req.user_id!);
      if (!projects) {
        console.error("getDeletedProjectsOfUsers returned null/undefined");
        projects = [];
      }
      return res.status(200).json(projects);
    } catch (err) {
      console.error(err);
      return res.status(500).send("Internal server error");
    }
  });
  return router;
};

const restore_deleted_project = (): Router => {
  let router = express.Router();
  router.post(
    "/deletedprojects/:id/restore",
    check_auth(),
    validationMdw,
    async (req: ExtRequest, res) => {
      try {
        let project = (
          await DeletedProjectService.getDeletedProjectById(req.params.id)
        )?.toObject();
        if (!project) {
          return res.status(404).send("Project not found");
        }
        let { deletedAt, ...restoredProject } = project;
        await ProjectService.addProject(restoredProject);
        await DeletedProjectService.deleteDeletedProjectById(project._id);
        return res.status(200).send("OK");
      } catch (err) {
        console.error(err);
        return res.status(500).send("Internal server error");
      }
    }
  );
  return router;
};

const cleanup_old_projects = (): Router => {
  let router = express.Router();
  router.get("/clean-recyclebin", async (req: ExtRequest, res) => {
    try {
      let oldprojects = await DeletedProjectService.getProjectsOverLifetime();
      for (const project of oldprojects) {
        try_delete_object_from_s3(result_path(project.id));
        try_delete_object_from_s3(result_model_path(project.id));
        try_delete_object_from_s3(query_path(project.id));
      }
      let ids = oldprojects.map((p) => p.id);
      await DeletedProjectService.deleteProjectsByIds(ids);
      return res.status(200).send("OK");
    } catch (err) {
      console.error(err);
      return res.status(500).send("Internal server error");
    }
  });
  return router;
};
function try_delete_object_from_s3(key: string) {
  let params: DeleteObjectRequest = {
    Bucket: process.env.S3_BUCKET_NAME,
    Key: key,
  };
  s3.deleteObject(params)
    .promise()
    .catch((e) => {
      console.error(`Error while deleting ${key}`);
      console.error(e);
    });
}

export {
  get_projects,
  get_userProjects,
  get_project_by_id,
  get_users_projects,
  update_project_results,
  delete_project,
  get_deleted_projects,
  restore_deleted_project,
  cleanup_old_projects,
};
