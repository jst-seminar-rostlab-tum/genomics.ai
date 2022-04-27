import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import express, { Response } from "express";
import processImageUpload from "../../../util/processImageUpload";
import InstitutionService from "../../../database/services/institution.service";
import { UpdateInstitutionDTO } from "../../../database/dtos/institution.dto";
import { IInstitution } from "../../../database/models/institution";

async function checkForValidValues(req: ExtRequest, res: Response): Promise<{ instId: string; bucket: string } | null> {
    let instId = req.params.id;
    let institution = instId ? await InstitutionService.getInstitutionById(instId) : null;
    if (!institution) {
        res.status(404).send("Institution not found");
        return null;
    }
    if (!institution.adminIds.includes(req.user_id!)) {
        res.status(403).send("Forbidden");
        return null;
    }
    if (!process.env.S3_PICTURE_BUCKET_NAME) {
        res.status(500).send("S3-BucketName is not set");
        return null;
    }
    return { instId, bucket: process.env.S3_PICTURE_BUCKET_NAME };
}

export function upload_institution_profilepicture_route() {
    let router = express.Router();
    router.post("/institutions/:id/profilepicture", check_auth(), async (req: ExtRequest, res) => {
        try {
            let valid = await checkForValidValues(req, res);
            if (!valid) return;
            let { instId, bucket } = valid;
            let result = await processImageUpload(req, { maxWidth: 400, maxHeight: 400, forceAspectRatio: true }, bucket, `inst_${instId}`);
            if (result.success === true) {
                const institutionUpdate: UpdateInstitutionDTO = {
                    profilePictureURL: result.objectUrl
                };
                await InstitutionService.updateInstitution(instId, institutionUpdate);
                res.status(200).contentType("text/plain").send(result.objectUrl);
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
            let valid = await checkForValidValues(req, res);
            if (!valid) return;
            let { instId, bucket } = valid;
            let result = await processImageUpload(req, { maxWidth: 1920, maxHeight: 1920, forceAspectRatio: false }, bucket, `instback_${instId}`);
            if (result.success === true) {
                const institutionUpdate: UpdateInstitutionDTO = {
                    backgroundPictureURL: result.objectUrl
                };
                await InstitutionService.updateInstitution(instId, institutionUpdate);
                res.status(200).contentType("text/plain").send(result.objectUrl);
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
