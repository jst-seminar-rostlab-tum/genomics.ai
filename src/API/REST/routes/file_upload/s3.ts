import AWS, {S3} from "aws-sdk";

const S3_OPTIONS: S3.Types.ClientConfiguration = {
    accessKeyId: process.env.S3_ACCESS_KEY_ID,
    secretAccessKey: process.env.S3_SECRET_ACCESS_KEY,
    endpoint: process.env.S3_ENDPOINT,
    s3ForcePathStyle: true,
    signatureVersion: 'v4'
};
const s3 = new AWS.S3(S3_OPTIONS);

export default s3;