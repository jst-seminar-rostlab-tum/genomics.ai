import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import { S3 } from "aws-sdk";
import s3 from "../../../util/s3";
import express from "express";
import UserService from "../../../database/services/user.service";
import { UpdateUserDTO } from "../../../database/dtos/user.dto";

export default function upload_user_avatar_route() {
    let router = express.Router();
    router.post("/upload_user_avatar", check_auth(), async (req: ExtRequest, res) => {
        let data = req.body;
        try {
            //TODO: check if data is valid image
            //TODO: handle different image files
            if (process.env.S3_PICTURE_BUCKET_NAME) {
                if (req.user_id) {
                    const params = {
                        Bucket: process.env.S3_PICTURE_BUCKET_NAME,
                        Key: `user_${req.user_id}.png`,
                        Body: data,
                        ContentType: "image/png"
                    };
                    let objectUrl = await new Promise<string>(function (resolve, reject) {
                        s3.upload(params, undefined, function (err, uploadData) {
                            if (err) {
                                console.error(err, err.stack || "Error while saving avatar");
                                reject(err);
                            } else {
                                resolve(uploadData.Location);
                            }
                        });
                    });
                    const userUpdate: UpdateUserDTO= {
                        avatarUrl: objectUrl
                    };
                    await UserService.updateUser(req.user_id,userUpdate)
                    res.status(200).send(objectUrl);
                } else {
                    res.status(401).send("Unauthorized");
                }
            } else res.status(500).send("S3-BucketName is not set");
        } catch (err) {
            console.log(err);
            res.status(500).send(err);
        }
    });
    return router;
}
