import e, { Response } from "express";
import { ExtRequest } from "../definitions/ext_request";
import sharp from "sharp";
import s3 from "./s3";
import {ImageUploadResult} from "../definitions/image_upload_result";

export default async function processImageUpload(
    req: ExtRequest,
    width: number,
    height: number,
    bucketName: string,
    resultFilename: string
): Promise<ImageUploadResult> {
    //If yes, the actual Content-Type is returned
    let contentType = req.is(["image/png", "image/jpeg"]);
    if (!contentType) {
        return {
            success: false,
            status: 400,
            message: "Invalid image format, only PNG and JPG are supported",
            error: new Error(`Invalid image format: ${req.headers["content-type"]}`)
        };
    }
    let croppedImage: Buffer, format: string;
    try {
        const image = sharp(req.body);
        const metadata = await image.metadata();
        const targetRatio = width / height;
        const maxWidth = width;
        const maxHeight = height;
        const sourceWidth = metadata.width!;
        const sourceHeight = metadata.height!;


        //Calculate target width and height, so the resulting picture is cropped to the requested aspect ratio
        //Afterwards, if the image is still to large it will be scaled down to the requested width and height
        //If the image is smaller, it will thus only be cropped to the requested aspect ratio.
        let targetWidth = targetRatio * sourceHeight;
        let targetHeight = sourceHeight;
        //Crop this amount left and right
        let cropX = (sourceWidth - targetWidth) / 2;
        let cropY = 0;

        if (sourceWidth < targetWidth) {
            targetWidth = sourceWidth;
            targetHeight = sourceWidth / targetRatio;
            cropX = 0;
            cropY = (sourceHeight - targetHeight) / 2;
        }

        let result = await image
            .extract({ left: Math.floor(cropX), top: Math.floor(cropY), width: Math.ceil(targetWidth), height: Math.ceil(targetHeight) })
            .resize(maxWidth, maxHeight, {
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
