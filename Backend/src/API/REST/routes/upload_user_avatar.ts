import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import { S3 } from "aws-sdk";
import s3 from "../../../util/s3";
import express from "express";
import UserService from "../../../database/services/user.service";
import { UpdateUserDTO } from "../../../database/dtos/user.dto";
import sharp from "sharp";
import processImageUpload from "../../../util/processImageUpload";

export default function upload_user_avatar_route() {
    let router = express.Router();
    router.post("/user-avatar", check_auth(), async (req: ExtRequest, res) => {
        try {
            if (process.env.S3_PICTURE_BUCKET_NAME) {
                    try {
                        let result = await processImageUpload(
                            req,
                            { maxWidth: 400, maxHeight: 400, forceAspectRatio: true },
                            process.env.S3_PICTURE_BUCKET_NAME,
                            `user_${req.user_id}`
                        );
                        if (result.success === true) {
                            const userUpdate: UpdateUserDTO = {
                                avatarUrl: result.objectUrl
                            };
                            await UserService.updateUser(req.user_id!, userUpdate);
                            res.status(200).contentType("text/plain").send(result.objectUrl);
                        } else {
                            const { status, message, error } = result;
                            res.status(status).send(message);
                            if (error) throw error;
                        }
                    } catch (e) {
                        console.error(e);
                    }
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    });
    return router;
}
