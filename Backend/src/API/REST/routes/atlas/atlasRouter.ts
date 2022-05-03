import express, { Router } from "express";
import AtlasService from "../../../../database/services/atlas.service";
import s3 from "../../../../util/s3";
import { validationMdw } from "../../middleware/validation";

/**
 *  Get details about an atlas.
 */
const get_atlas = (): Router => {
  let router = express.Router();

  router.get("/atlas/:id", validationMdw, async (req: any, res) => {
    const atlasId: string = req.params.id;

    try {
      const atlas = await AtlasService.getAtlasById(atlasId);
      return res.status(200).json(atlas);
    } catch (err) {
      console.error("Error getting information about the atlas!");
      console.error(JSON.stringify(err));
      console.error(err);
      return res.status(500).send("Unable to retrieve information about the atlas.");
    }
  });
  return router;
};

const get_atlas_visualization = (): Router => {
  let router = express.Router();

  router.get("/atlas/:id/visualization", validationMdw, async (req: any, res) => {
    //TODO: Using presigned urls at the moment, instead of a public bucket, is a temporary solution for the moment.
    const atlasId = req.params.id;

    try {
      const atlas = await AtlasService.getAtlasById(atlasId);
      if (!atlas) return res.status(404).send("Atlas not found");
      let params: any = {
        Bucket: process.env.S3_BUCKET_NAME!,
        Key: `atlas/${atlasId}/visualization.csv`,
        Expires: 60 * 60 * 24 * 7 - 1, // one week minus one second
      };
      let presignedUrl = await s3.getSignedUrlPromise("getObject", params);
      return res.status(200).contentType("text/plain").send(presignedUrl);
    } catch (err) {
      console.error(err);
      return res.status(500).send("Internal error");
    }
  });
  return router;
};

/**
 *  Get all available atlases.
 */
const get_allAtlases = (): Router => {
  let router = express.Router();

  router.get("/atlases", validationMdw, async (req: any, res) => {
    try {
      const atlases = await AtlasService.getAllAtlases();
      return res.status(200).json(atlases);
    } catch (err) {
      console.error("Error accessing the atlases!");
      console.error(JSON.stringify(err));
      console.error(err);
      return res.status(500).send("Unable to access the atlases.");
    }
  });
  return router;
};

export { get_atlas, get_atlas_visualization, get_allAtlases };
