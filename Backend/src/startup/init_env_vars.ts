import dotenv from 'dotenv';

dotenv.config();

export function init_env_vars() {
    // check for required env-vars
    const required_env_vars = [
        "DATABASE_URI",
        "CLOUD_RUN_URL",
        "S3_ACCESS_KEY_ID",
        "S3_SECRET_ACCESS_KEY",
        "S3_ENDPOINT",
        "S3_BUCKET_NAME",
        "MAILGUN_API_KEY",
        "MAIL_DOMAIN"
    ];
    required_env_vars.forEach((required_env_var) => {
        if (!(required_env_var in process.env))
            console.warn("WARNING: the environment variable '" + required_env_var + "' is not defined!");
    })

    // set default values for missing env vars
    function setStdEnvValue(env: string, value: any) {
        if (!process.env.hasOwnProperty(env))
            process.env[env] = value;
    }

    setStdEnvValue("PORT", "8050");
    setStdEnvValue("DATABASE_URI", "mongodb://localhost:27017/dev");

    setStdEnvValue("S3_ENDPOINT", "http://127.0.0.1:9000");
    setStdEnvValue("S3_BUCKET_NAME", "vaultgovsg");
}