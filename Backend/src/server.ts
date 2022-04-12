//Setup environment first
import * as Startup from "./startup/_startup";
Startup.init_environment();
Startup.init_env_vars();

import {HTTP_Server} from "./http_server/http_server";
import API_Host from "./API/_api_host";
import {Database} from "./database/database";

console.log(" *** Startup *** ");

let http_server = new HTTP_Server();
let apiHost = new API_Host(http_server);
let database = new Database();

(async ()=>{
    // connect to database
    await database.connect();
    // initialise APIs
    await http_server.setup();
    await apiHost.init();
})().then(()=>{
    console.log(" *** Startup complete. *** ");
})
