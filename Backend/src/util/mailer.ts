import path from "path";

import nodemailer from "nodemailer";
import mailgunTransport from "nodemailer-mailgun-transport";
import fs from "fs/promises";
import handlebars from "handlebars";

class Mailer {
    transport: nodemailer.Transporter | null;

    constructor() {
        if (!process.env.MAILGUN_API_KEY || process.env.MAILGUN_API_KEY == "") {
            this.transport = null;
        } else {
            this.transport = nodemailer.createTransport(
                mailgunTransport({
                    auth: {
                        api_key: process.env.MAILGUN_API_KEY!,
                        domain: process.env.MAIL_DOMAIN
                    },
                    host: process.env.MAILGUN_HOST
                })
            );
        }
    }

    async send(to : string, subject : string, template_name : string, data : any) {
        let self = this;
        let texttemplate = await fs.readFile(path.join(__dirname, "./../../src/views/mails", template_name, "text.txt"), 'utf-8')
        let htmltemplate = await fs.readFile(path.join(__dirname, "./../../src/views/mails", template_name, "html.html"), 'utf-8');

        let rendered_txt = handlebars.compile(texttemplate)(data);
        let rendered_html = handlebars.compile(htmltemplate)(data);

        let mail = {
            from: `GeneCruncher <info@${process.env.MAIL_DOMAIN}>`,
            to: to,
            subject: subject,
            text: rendered_txt,
            html: rendered_html
        };
        if (!self.transport) {
            console.log("No mailgun key provided: logging email here:");
            console.log(`From: ${mail.from}`);
            console.log(`To: ${mail.to}`);
            console.log(`Subject: ${mail.subject}`);
            console.log(`Text: ${mail.text}`);
            console.log(`HTML: ${mail.html}`);
        } else {
            return self.transport.sendMail(mail);
        }
    }

    async send_verification_mail(firstname: string, recipient: string, token: string) {
        // return this.send(recipient, "Verify Email", "signup_confirm_email", {link: `https://api.genecruncher.com/verify/${token}`, firstname:firstname})
        return this.send(recipient, "Verify Email", "signup_confirm_email", {link: `${process.env.API_URL}/verify/${token}`, firstname:firstname})
    }
}

export const mailer = new Mailer();
