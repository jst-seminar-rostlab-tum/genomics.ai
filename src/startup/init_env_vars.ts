export function init_env_vars (){
    // check for required env-vars
    const required_env_vars = [
        "test_env_var"
    ];
    required_env_vars.forEach((required_env_var)=>{
        if(!(required_env_var in process.env))
            console.warn("WARNING: the environment variable '" + required_env_var + "' is not defined!");
    })

    // set default values for missing env vars
    function setStdEnvValue(env:string, value:any){
        if(!process.env.hasOwnProperty(env))
            process.env[env] = value;
    }

    setStdEnvValue("PORT", "8050");
}