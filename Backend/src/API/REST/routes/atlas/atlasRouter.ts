import express, {Router} from "express";
import AtlasService from "../../../../database/services/atlas.service";

/**
 *  Get details about an atlas.
 */
const get_atlas = () : Router => {
    let router = express.Router();

    router.get("/atlas/:id", async (req: any, res) => {
        const atlasId: string = req.params.id;
        if (!atlasId)
            return res.status(400).send("Missing atlas id");

        try {
            const atlas = AtlasService.getAtlasById(atlasId);
            return res.status(200).json(atlas);
        } catch(err) {
            console.error("Error getting information about the atlas!");
            console.error(JSON.stringify(err));
            console.error(err);
            return res.status(500).send("Unable to retrieve information about the atlas.");
        }
    })
    return router;
}

export { get_atlas };
