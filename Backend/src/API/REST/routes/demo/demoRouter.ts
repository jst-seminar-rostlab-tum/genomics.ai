import express, { Router } from 'express';
import DemoService from "../../../../database/services/demo.service";

/**
 * Return details about all the demo datasets that are available
 */
const get_allDemos = (): Router => {
    let router = express.Router();

    router.get("/demos", async (req: any, res) => {
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