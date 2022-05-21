import express, { Router } from "express";
import ModelService from "../../../../database/services/model.service";

/**
 *  Get details about a model.
 */
const get_model = (): Router => {
  let router = express.Router();

  router.get("/model/:id", async (req: any, res) => {
    const modelId: string = req.params.id;
    if (!modelId) return res.status(400).send("Missing model id");

    try {
      const model = await ModelService.getModelById(modelId);
      return res.status(200).json(model);
    } catch (err) {
      console.error("Error getting information about the model!");
      console.error(JSON.stringify(err));
      console.error(err);
      return res.status(500).send("Unable to retrieve information about the model.");
    }
  });
  return router;
};

/**
 *  Get all available models.
 */
const get_allModels = (): Router => {
  let router = express.Router();

  router.get("/models", async (req: any, res) => {
    try {
      const models = await ModelService.getAllModels();
      return res.status(200).json(models);
    } catch (err) {
      console.error("Error accessing the models!");
      console.error(JSON.stringify(err));
      console.error(err);
      return res.status(500).send("Unable to access the models.");
    }
  });
  return router;
};

export { get_model, get_allModels };
