import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import express, { Response } from "express";
import processImageUpload from "../../../util/processImageUpload";
import InstitutionService from "../../../database/services/institution.service";
import { UpdateInstitutionDTO } from "../../../database/dtos/institution.dto";
import { IInstitution } from "../../../database/models/institution";
import s3 from "../../../util/s3";
import { S3 } from "aws-sdk";

async function retrieveInstitutionFromRequest(
  req: ExtRequest,
  res: Response
): Promise<string | null> {
  let instId = req.params.id;
  let institution = instId ? await InstitutionService.getInstitutionById(instId) : null;
  if (!institution) {
    res.status(404).send("Institution not found");
    return null;
  }
  if (!institution.adminIds.includes(req.user_id!)) {
    res.status(401).send("Unauthorized");
    return null;
  }
  return instId;
}

async function reset_picture_logic(
  req: ExtRequest,
  res: Response,
  updateDatabase: (instId: string) => Promise<string | null | undefined>
): Promise<any> {
  try {
    let instId = await retrieveInstitutionFromRequest(req, res);
    if (!instId) return;

    if (!process.env.S3_PICTURE_BUCKET_NAME) {
      return res.status(500).send("S3-BucketName is not set");
    }

    let oldUrl = await updateDatabase(instId);
    if (oldUrl) {
      let key = oldUrl.substring(oldUrl.lastIndexOf("/") + 1, oldUrl.lastIndexOf("?"));
      let params: S3.Types.DeleteObjectRequest = {
        Bucket: process.env.S3_PICTURE_BUCKET_NAME,
        Key: key,
      };
      await s3.deleteObject(params).promise();
    }
    res.status(200).send("OK");
  } catch (err) {
    console.log(err);
    res.status(500).send(err);
  }
}

export function reset_institution_profilepicture_route() {
  let router = express.Router();
  router.delete("/institutions/:id/profilepicture", check_auth(), async (req: ExtRequest, res) => {
    await reset_picture_logic(req, res, async (instId) => {
      return await InstitutionService.unsetProfilePicture(instId);
    });
  });
  return router;
}

export function reset_institution_backgroundpicture_route() {
  let router = express.Router();
  router.delete(
    "/institutions/:id/backgroundpicture",
    check_auth(),
    async (req: ExtRequest, res) => {
      await reset_picture_logic(req, res, async (instId) => {
        return await InstitutionService.unsetBackgroundPicture(instId);
      });
    }
  );
  return router;
}
