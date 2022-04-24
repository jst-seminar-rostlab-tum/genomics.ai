import e, { Response } from "express";
import { ExtRequest } from "../definitions/ext_request";
import sharp from "sharp";
import s3 from "./s3";

export default async function processImageUpload(
    req: ExtRequest,
    width: number,
    height: number,
    bucketName: string,
    resultFilename: string
): Promise<ImageUploadResult> {
    let contentType = req.is(["image/png", "image/jpeg"]);
    if (!contentType) {
        return { success: false, status: 400, message: "Invalid image format, only PNG and JPG are supported", error: new Error(`Invalid image format: ${req.headers["content-type"]}`) };
    }
    let croppedImage: Buffer, format: string;
    try {
        const result = await sharp(req.body)
            .resize(width, height, {
                fit: "cover",
                withoutEnlargement: true
            })
            .toBuffer({ resolveWithObject: true });
        
        croppedImage = result.data;
        format = result.info.format;
    } catch (e) {
        return { success: false, status: 500, message: "Image processing failed", error: e };
    }

    let extension;
    if (format == "png") {
        extension = ".png";
    } else if (format == "jpeg") {
        extension = ".jpg";
    } else {
        return { success: false, status: 500, message: "Image processing failed" };
    }

    try {
        const params = {
            Bucket: bucketName,
            Key: `${resultFilename}${extension}`,
            Body: croppedImage,
            ContentType: contentType
        };
        let objectUrl = await new Promise<string>(function (resolve, reject) {
            s3.upload(params, undefined, function (err, uploadData) {
                if (err) {
                    reject(err);
                } else {
                    resolve(uploadData.Location);
                }
            });
        });
        return { success: true, objectUrl };
    } catch (e) {
        return { success: false, status: 500, message: "Error while saving image", error: e };
    }
}
