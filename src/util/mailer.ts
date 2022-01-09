import path from "path";

import nodemailer from "nodemailer";
import mailgunTransport from "nodemailer-mailgun-transport";
import fs from "fs/promises";
import handlebars from "handlebars";

class Mailer {
    transport : nodemailer.Transporter;

    constructor(){
        this.transport = nodemailer.createTransport(mailgunTransport({
            auth: {
                api_key: process.env.MAILGUN_API_KEY!,
                domain: process.env.MAIL_DOMAIN
            }
        }));
    }

    async send(to : string, subject : string, template_name : string, data : any) {
        let self = this;
        let texttemplate = await fs.readFile(path.join(__dirname, "./../../views/mails", template_name, "text.txt"), 'utf-8')
        let htmltemplate = await fs.readFile(path.join(__dirname, "./../../views/mails", template_name, "html.html"), 'utf-8');

        let rendered_txt = handlebars.compile(texttemplate)(data);
        let rendered_html = handlebars.compile(htmltemplate)(data);
        return self.transport.sendMail({
            from: `GeneCruncher <info@${process.env.MAIL_DOMAIN}>`,
            to: to,
            subject: subject,
            text: rendered_txt,
            html: rendered_html
        });
    }

    async send_verification_mail(firstname: string, recipient: string, token: string) {
        return this.send(recipient, "Verify Email", "signup_confirm_email", {link: `http://localhost:8050/verify/${token}`, firstname:firstname})
    }
}

export const mailer = new Mailer();