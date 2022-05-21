import check_auth from "../middleware/check_auth";
import { ExtRequest } from "../../../definitions/ext_request";
import express from "express";
import UserService from "../../../database/services/user.service";
import { userModel } from "../../../database/models/user";
import s3 from "../../../util/s3";
import { S3 } from "aws-sdk";

export default function reset_user_avatar_route() {
  let router = express.Router();
  router.delete("/user-avatar", check_auth(), async (req: ExtRequest, res) => {
    try {
      if (!process.env.S3_PICTURE_BUCKET_NAME) {
        res.status(500).send("S3-BucketName is not set");
        return null;
      }
      let oldUrl = await UserService.unsetUserAvatar(req.user_id!);

      if (oldUrl) {
        let key = oldUrl.substring(oldUrl.lastIndexOf("/") + 1, oldUrl.lastIndexOf("?"));
        console.log(key);
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
  });
  return router;
}
