import express, { Router } from 'express';
import check_auth from '../../middleware/check_auth';
//import DemoService from "../../../../database/demo.service";
// TODO: write the code for the demo service

/**
 * Return details about all the demo datasets that are available
 */
const get_allDemos = (): Router => {
    let router = express.Router();

    router.get("/demos", check_auth, async (req: any, res) => {
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

// todo: the demo dataset object should look like this

// The demo dataset object should look like: 
// The demo dataset object should be: 
// {
// demoId
// title
// atlas
// model
// url to fetch the csv from 
// }

// TODO: Create the model for mongoDB, i.e. the schema and teh service. 

// TODO: Understand how the csv are stored and fetched. 

// TODO: set up the server locally, or at least only the s3 bucket I believe. Necessary in order to test the csv's 