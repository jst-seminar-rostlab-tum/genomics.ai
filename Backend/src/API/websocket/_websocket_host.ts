import WebSocket from "ws";
import { setup_ws_server } from "./ws_setup/_setup_ws_server";
import { ws_send } from "./ws_send";
import { HTTP_Server } from "../../http_server/http_server";

export default class Websocket_Host {
  // ---- fields ----
  // instance of the websocket-server running through the http_server
  protected wss: WebSocket.Server | null = null;

  // instance of the http-server hosting the websocket connections
  protected http_server: HTTP_Server;

  // interval for periodically pinging all connections (instantiated through setInterval in setup_ws_ping_interval)
  protected ws_ping_interval: NodeJS.Timeout | null = null;
  // ------------

  // ---- methods ----
  // send a message to a connected websocket
  protected ws_send = ws_send;
  // ------------

  constructor(http_server: HTTP_Server) {
    this.http_server = http_server;
  }

  // launches the websocket-server and performs all necessary setup
  public async init() {
    await setup_ws_server.call(this);

    console.log(`Websocket-Server operational.`);
    return this;
  }
}
