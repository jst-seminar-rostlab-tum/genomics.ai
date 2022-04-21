import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import { S3 } from "aws-sdk";
import s3 from "../../../util/s3";
import express from "express";
import { userModel } from "../../../database/models/user";

export default function upload_user_avatar_route() {
    let router = express.Router();
    router.post("/upload_user_avatar", check_auth(), async (req: ExtRequest, res) => {
        let data = req.body;
        console.log(data);
        try {
            if (process.env.S3_PICTURE_BUCKET_NAME) {
                if (req.user_id) {
                    const params = {
                        Bucket: process.env.S3_PICTURE_BUCKET_NAME,
                        Key: `user_${req.user_id}.png`,
                        Body: data,
                        ContentType: "image/png"
                    };
                    let objectUrl = await new Promise(function (resolve, reject) {
                        s3.upload(params, undefined, function (err, uploadData) {
                            if (err) {
                                console.error(err, err.stack || "Error while saving avatar");
                                reject(err);
                            } else {
                                resolve(uploadData.Location);
                            }
                        });
                    });
                    await userModel.updateOne({_id: req.user_id}, {
                        $set: {
                            avatarUrl: objectUrl
                        }
                    })
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
