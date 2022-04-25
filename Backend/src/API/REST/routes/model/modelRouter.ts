import express, {Router} from "express";
import {ExtRequest} from "../../../../definitions/ext_request";
import ModelService from "../../../../database/services/model.service";

/**
 *  Get details about a model.
 */
const get_model_details = () : Router => {
    let router = express.Router();

    router.get("/model/:id", async (req: any, res) => {
        const modelId: string = req.params.id;
        if (!modelId)
            return res.status(400).send("Missing model id");

        try {
            const model = ModelService.getModelById(modelId);
            return res.status(200).json(model);
        } catch(err) {
            console.error("Error getting information about the model!");
            console.error(JSON.stringify(err));
            console.error(err);
            return res.status(500).send("Unable to retrieve information about the model.");
        }
    })
    return router;
}

export { get_model_details };
