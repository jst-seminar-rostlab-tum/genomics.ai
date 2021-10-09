import {HTTP_Server} from "./http_server/http_server";
import * as Startup from "./startup/_startup";
import API_Host from "./API/_api_host";

console.log(" *** Startup *** ");
Startup.init_environment();
Startup.init_env_vars();

let http_server = new HTTP_Server();
let apiHost = new API_Host(http_server);

(async ()=>{
    // initialise APIs
    await http_server.setup();
    await apiHost.init();
})().then(()=>{
    console.log(" *** Startup complete. *** ");
})
