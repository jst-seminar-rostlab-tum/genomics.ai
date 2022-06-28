import express, { Router } from 'express';
import check_auth from '../../middleware/check_auth';
import DemoService from "../../../../database/services/demo.service";

/**
 * Return details about all the demo datasets that are available
 */
const get_allDemos = (): Router => {
    let router = express.Router();

    router.get("/demos", check_auth(), async (req: any, res) => {
        try {
            const models = await DemoService.getAllDemos();
            return res.status(200).json(models);
        } catch (err) {
            console.error("Error accesssing the demos");
            console.error(err);
            return res.status(500).send("Unable to access the demos.");
        }
    });
    return router;
}

export { get_allDemos };

// TODO: Understand how the csv's are stored and fetched. 

// TODO: set up the server locally, or at least only the s3 bucket I believe. Necessary in order to test the csv's 

// try out the demo projects locally, create th