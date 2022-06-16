import { Request } from "express";
import { ObjectId } from "mongoose";
import InstitutionService from "../../../database/services/institution.service";

export interface ExtInstRequest extends Request {
  user_id?: ObjectId; // declare optional property "user_id"
}

export const institution_admin_auth = (req: ExtInstRequest, res: any, next: any) => {
  const current_user = req.user_id;
  const institution_id = req.params.id;

  try {
    InstitutionService.getInstitutionById(institution_id)
      .then((inst) => {
        if (!inst) {
          return res.status(404).send("There is no institution with id: " + institution_id);
        }
        if (!inst.adminIds.includes(current_user!)) {
          return res.status(401).send("Invalid permissions to do this operations!");
        }
        next();
      })
      .catch((error) => {
        return res
          .status(500)
          .send("Error during authentication: Failed to fetch the institution, " + error);
      });
  } catch (e) {
    return res
      .status(500)
      .send("Error during authentication: Failed to check auth for institution");
  }
};

export const institution_member_auth = (req: ExtInstRequest, res: any, next: any) => {
  const current_user = req.user_id;
  const institution_id = req.params.id;

  try {
    InstitutionService.getInstitutionById(institution_id)
      .then((inst) => {
        if (!inst) {
          return res.status(404).send("There is no institution with id: " + institution_id);
        }
        if (!inst.memberIds.includes(current_user!)) {
          return res.status(401).send("Invalid permissions to do this operations!");
        }
        next();
      })
      .catch((error) => {
        return res
          .status(500)
          .send("Error during authentication: Failed to fetch the institution, " + error);
      });
  } catch (e) {
    return res
      .status(500)
      .send("Error during authentication: Failed to check auth for institution");
  }
};

export const institution_admin_or_member_auth = (req: ExtInstRequest, res: any, next: any) => {
  const current_user = req.user_id;
  const institution_id = req.params.id;

  try {
    InstitutionService.getInstitutionById(institution_id)
      .then((inst) => {
        if (!inst) {
          return res.status(404).send("There is no institution with id: " + institution_id);
        }
        if (!(inst.adminIds.includes(current_user!) || inst.memberIds.includes(current_user!))) {
          return res.status(401).send("Invalid permissions to do this operations!");
        }
        next();
      })
      .catch((error) => {
        return res
          .status(500)
          .send("Error during authentication: Failed to fetch the institution, " + error);
      });
  } catch (e) {
    return res
      .status(500)
      .send("Error during authentication: Failed to check auth for institution");
  }
};
