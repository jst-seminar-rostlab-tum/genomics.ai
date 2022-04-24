import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import express from "express";
import processImageUpload from "../../../util/processImageUpload";
import InstitutionService from "../../../database/services/institution.service";
import { UpdateInstitutionDTO } from "../../../database/dtos/institution.dto";
import {IInstitution} from "../../../database/models/institution";

async function tryGetInstitution(req: ExtRequest, res):Promise<{instId:string, institution:IInstitution}> {
    if (!req.user_id) {
        res.status(401).send("Not authenticated");
        return null;
    }
    let instId = req.params.id;
    let institution = instId ? await InstitutionService.getInstitutionById(instId) : null;
    if (!institution) {
        res.status(404).send("Institution not found");
        return null;
    }
    if (!institution.adminIds.includes(req.user_id)) {
        res.status(401).send("Unauthorized");
        return null;
    }
    if (!process.env.S3_PICTURE_BUCKET_NAME) {
        res.status(500).send("S3-BucketName is not set");
        return null;
    }
    return {instId,institution};
}

export function upload_institution_profilepicture_route() {
    let router = express.Router();
    router.post("/institutions/:id/profilepicture", check_auth(), async (req: ExtRequest, res) => {
        try {
            let {instId, institution} = await tryGetInstitution(req,res);
            let result = await processImageUpload(req, 400, 400, process.env.S3_PICTURE_BUCKET_NAME, `inst_${instId}`);
            if (result.success === true) {
                const institutionUpdate: UpdateInstitutionDTO = {
                    profilePictureURL: result.objectUrl
                };
                await InstitutionService.updateInstitution(req.user_id, institutionUpdate);
                res.status(200).send(result.objectUrl);
            } else {
                const { status, message, error } = result;
                if (error) console.log(error);
                return res.status(status).send(message);
            }
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    });
    return router;
}

export function upload_institution_backgroundpicture_route() {
    let router = express.Router();
    router.post("/institutions/:id/backgroundpicture", check_auth(), async (req: ExtRequest, res) => {
        try {
            let {instId, institution} = await tryGetInstitution(req,res);
            //TODO: change width+height to value discussed with frontend
            let result = await processImageUpload(req, 1600, 600, process.env.S3_PICTURE_BUCKET_NAME, `instback_${instId}`);
            if (result.success === true) {
                const institutionUpdate: UpdateInstitutionDTO = {
                    backgroundPictureUrl: result.objectUrl
                };
                await InstitutionService.updateInstitution(req.user_id, institutionUpdate);
                res.status(200).send(result.objectUrl);
            } else {
                const { status, message, error } = result;
                if (error) console.log(error);
                return res.status(status).send(message);
            }
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    });
    return router;
}
