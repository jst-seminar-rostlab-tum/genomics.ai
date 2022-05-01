import express, { Router } from "express";
import AtlasService from "../../../../database/services/atlas.service";

/**
 *  Get details about an atlas.
 */
const get_atlas = (): Router => {
  let router = express.Router();

  router.get("/atlas/:id", async (req: any, res) => {
    const atlasId: string = req.params.id;
    if (!atlasId) return res.status(400).send("Missing atlas id");

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

/**
 *  Get all available atlases.
 */
const get_allAtlases = (): Router => {
  let router = express.Router();

  router.get("/atlases", async (req: any, res) => {
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

export { get_atlas, get_allAtlases };
