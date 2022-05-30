import AWS, { S3 } from "aws-sdk";
import { DeleteObjectRequest } from "aws-sdk/clients/s3";

const S3_OPTIONS: S3.Types.ClientConfiguration = {
  accessKeyId: process.env.S3_ACCESS_KEY_ID,
  secretAccessKey: process.env.S3_SECRET_ACCESS_KEY,
  endpoint: process.env.S3_ENDPOINT,
  s3ForcePathStyle: true,
  signatureVersion: "v4",
};
const s3 = new AWS.S3(S3_OPTIONS);

export function try_delete_object_from_s3(key: string) {
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

export default s3;
