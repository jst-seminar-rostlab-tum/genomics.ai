import e, { Response } from "express";
import { ExtRequest } from "../definitions/ext_request";
import Jimp from "jimp";
import s3 from "./s3";
import { ImageUploadResult } from "../definitions/image_upload_result";

export default async function processImageUpload(
  req: ExtRequest,
  scalingOptions: { maxWidth: number; maxHeight: number; forceAspectRatio: boolean },
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
      error: new Error(`Invalid image format: ${req.headers["content-type"]}`),
    };
  }
  let croppedImage: Buffer;
  try {
    const image = await Jimp.read(req.body);
    const sourceWidth = image.getWidth();
    const sourceHeight = image.getHeight();
    let maxWidth = scalingOptions.maxWidth;
    let maxHeight = scalingOptions.maxHeight;
    let cropWidth = sourceWidth;
    let cropHeight = sourceHeight;
    let cropX = 0;
    let cropY = 0;
    if (scalingOptions.forceAspectRatio) {
      const targetRatio = maxWidth / maxHeight;

      //Calculate target width and height, so the resulting picture is cropped to the requested aspect ratio
      //Afterwards, if the image is still to large it will be scaled down to the requested width and height
      //If the image is smaller, it will thus only be cropped to the requested aspect ratio.
      cropWidth = targetRatio * sourceHeight;
      cropHeight = sourceHeight;
      //Crop this amount left and right
      cropX = (sourceWidth - cropWidth) / 2;
      cropY = 0;

      if (sourceWidth < cropWidth) {
        cropWidth = sourceWidth;
        cropHeight = sourceWidth / targetRatio;
        cropX = 0;
        cropY = (sourceHeight - cropHeight) / 2;
      }
    }

    let ops = image.crop(
      Math.floor(cropX),
      Math.floor(cropY),
      Math.ceil(cropWidth),
      Math.ceil(cropHeight)
    );
    if (cropWidth > maxWidth || cropHeight > maxHeight) {
      ops = ops.scaleToFit(maxWidth, maxHeight);
    }
    let buffer = await ops.getBufferAsync(contentType);

    croppedImage = buffer;
  } catch (e) {
    console.log(e);
    return { success: false, status: 500, message: "Image processing failed", error: e };
  }

  let extension;
  if (contentType == "image/png") {
    extension = ".png";
  } else if (contentType == "image/jpeg") {
    extension = ".jpg";
  } else {
    console.error("Wrong content type? " + contentType);
    return { success: false, status: 500, message: "Image processing failed" };
  }

  try {
    const params = {
      Bucket: bucketName,
      Key: `${resultFilename}`,
      Body: croppedImage,
      ContentType: contentType,
    };
    let uploadData = await s3.upload(params).promise();
    let objectUrl = `${uploadData.Location}?c=${Date.now()}`;
    return { success: true, objectUrl };
  } catch (e) {
    return { success: false, status: 500, message: "Error while saving image", error: e };
  }
}
