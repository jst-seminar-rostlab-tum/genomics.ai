import express, { Router } from "express";
import check_auth from "../../middleware/check_auth";
import {ExtRequest} from "../../../../definitions/ext_request";
import {mailer} from "../../../../util/mailer";

const contact_us = ():Router =>{
    let router = express.Router();
    router.post("/contact",  async (req: ExtRequest, res: any) => {
        try {
            const { email, firstName, lastName, message, subject } = req.body;

            await mailer.send(
                process.env.CONTACT_US,
                subject,
                "contact_us",
                {
                    email: email,
                    name: firstName + " " + lastName,
                    message: message,
                }
            );
            return res.status(200).send("Contact info sent successfully");
        } catch (e) {
            console.error(JSON.stringify(e));
            console.error(e);
            return res.status(500).send("Contact info could not be sent");
        }
    });
    return router;
}

export {contact_us}