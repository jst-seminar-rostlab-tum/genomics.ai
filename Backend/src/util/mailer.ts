import path from "path";

import nodemailer from "nodemailer";
import mailgunTransport from "nodemailer-mailgun-transport";
import fs from "fs/promises";
import handlebars from "handlebars";
import { v4 as uuid } from "uuid";

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
            domain: process.env.MAIL_DOMAIN,
          },
          host: process.env.MAILGUN_HOST,
        })
      );
    }
  }

  async send(
    to: string,
    subject: string,
    template_name: string,
    data: any,
    attachments?: { [key: string]: string }
  ) {
    let self = this;
    let texttemplate = await fs.readFile(
      path.join(__dirname, "./../views/mails", template_name, "text.txt"),
      "utf-8"
    );
    let htmltemplate = await fs.readFile(
      path.join(__dirname, "./../views/mails", template_name, "index.html"),
      "utf-8"
    );

    let attachmentsWithCid = [];
    if (attachments) {
      for (const [ident, filepath] of Object.entries(attachments)) {
        let id = uuid();
        let filename = path.basename(filepath);
        let cid = `${id}.${filename}`;
        data[ident] = `cid:${cid}`;

        let fullpath = path.join(__dirname, "./../views/mails", template_name, "front", filepath);
        let content = await fs.readFile(fullpath);

        //mailgun uses cid as filename
        attachmentsWithCid.push({
          cid: cid,
          content: content,
        });
      }
    }

    let rendered_txt = handlebars.compile(texttemplate)(data);
    let rendered_html = handlebars.compile(htmltemplate)(data);

    console.log(attachmentsWithCid);

    let mail = {
      from: `GeneCruncher <noreply@${process.env.MAIL_DOMAIN}>`,
      to: to,
      subject: subject,
      text: rendered_txt,
      html: rendered_html,
      attachments: attachmentsWithCid,
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
    return this.send(
      recipient,
      "Verify Email",
      "signup_confirm_email",
      {
        link: `${process.env.API_URL}/verify/${token}`,
        firstname: firstname,
        genecruncher_logo: `${process.env.STATIC_FILES_URL}/genecruncher_logo_email.svg`,
      },
    );
  }

  async send_password_reset_request_mail(firstname: string, recipient: string, token: string) {
    return this.send(
      recipient,
      "[GeneCruncher] Please reset your password",
      "password_reset_request_email",
      {
        firstname: firstname,
        link: `${process.env.FRONTEND_URL}/password_reset?token=${token}`,
        new_reset_link: `${process.env.FRONTEND_URL}/password_reset`,
        genecruncher_logo: `${process.env.STATIC_FILES_URL}/genecruncher_logo_email.svg`,
      },
    );
  }
  async send_password_reset_confirmation_mail(firstname: string, recipient: string) {
    return this.send(
      recipient,
      "[GeneCruncher] Your password was reset successfully",
      "password_reset_confirmation_email",
      {
        firstname: firstname,
        email: recipient,
        link: `${process.env.FRONTEND_URL}/password_reset`,
        genecruncher_logo: `${process.env.STATIC_FILES_URL}/genecruncher_logo_email.svg`,
      },
    );
  }
  async send_invitation_to_institution_mail(
    firstname: string,
    recipient: string,
    institution: string,
    country: string,
    link: string
  ) {
    return this.send(
      recipient,
      "[GeneCruncher] Invitation to an institution",
      "invitation_to_institution",
      {
        institution,
        country,
        firstname,
        link,
        genecruncher_logo: `${process.env.STATIC_FILES_URL}/genecruncher_logo_email.svg`,
      },
    );
  }
  async send_invitation_to_team_mail(
    firstname: string,
    recipient: string,
    teamname: string,
    link: string
  ) {
    return this.send(
      recipient,
      "[GeneCruncher] Invitation to a team",
      "invitation_to_team",
      {
        firstname,
        teamname,
        link,
        genecruncher_logo: `${process.env.STATIC_FILES_URL}/genecruncher_logo_email.svg`,
      }
    );
  }
}

export const mailer = new Mailer();
