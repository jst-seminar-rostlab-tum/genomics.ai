import Websocket_Host from "./websocket/_websocket_host";
import {HTTP_Server} from "../http_server/http_server";
import REST_Host from "./REST/_rest_host";

export default class API_Host {
    // constructor-provided references to the http_server to be used by this API_Host
    protected http_server : HTTP_Server;

    // child-classes
    protected websocket_host : Websocket_Host;
    protected rest_host : REST_Host;

    constructor(http_server : HTTP_Server) {
        this.http_server = http_server;
        this.websocket_host = new Websocket_Host(http_server);
        this.rest_host = new REST_Host(http_server);
    }

    public async init () {
        await this.websocket_host.init();
        await this.rest_host.init();

        return this;
    }
}