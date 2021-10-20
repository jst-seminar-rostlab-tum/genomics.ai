import REST_Host from "./_rest_host";
import express, {Router} from "express";
import rateLimit from "express-rate-limit";
import bodyParser from "body-parser";

// express.js middleware before requests hit route-handlers
export function express_middleware(this:REST_Host) : Router {
    let router = express.Router();

    router.use(rateLimit({
        windowMs: 60 * 1000, // 1min
        max: 1000, // limit each IP to 1000 requests per windowMs
        message: "Limit exceeded. Try again in 3 hours."
    }))

    router.use(bodyParser.json());
    router.use(bodyParser.urlencoded({ extended: true }));

    return router;
}