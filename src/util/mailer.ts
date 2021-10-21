

export function send_mail(recipient: string, subject: string, body: string) {
    const API_KEY = process.env.MAIL_API_KEY;
    const DOMAIN = process.env.MAIL_DOMAIN;
    /*const mailer = new NodeMailgun(API_KEY, DOMAIN);

    mailer.fromEmail = "noreply@genecruncher.io";
    mailer.fromTitle = subject;

    mailer
        .init()
        .send(recipient, subject, body)
        .then((_:any) => console.log(`Mail sent to ${recipient}`))
        .catch((error) => console.error(error));*/
}

export function send_verification_mail(recipient: string, token: string) {
    send_mail(
        recipient, 
        "Verification",
        `Thank you for signing up! Please verify your e-mail address using this link: http://localhost:8050/verify/${token}.\nThank you!`
    );
}