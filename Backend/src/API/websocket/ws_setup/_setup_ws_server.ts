import Websocket_Host from "../_websocket_host";
import WebSocket from "ws";
import { ExtWebSocket } from "../../../definitions/ext_web_socket";
import { on_message } from "../ws_on_message/_ws_on_message";
import { on_close } from "../ws_on_close";
import { setup_ws_ping_interval } from "./setup_ws_ping_interval";

// setup the websocket-server on top of the http_server
export function setup_ws_server(this: Websocket_Host) {
  const server = this.http_server.server;

  this.wss = new WebSocket.Server({ server });
  this.wss?.on("connection", (ws: ExtWebSocket) => {
    ws.on("pong", () => (ws.timeout = false));
    ws.on("close", () => on_close.call(this, ws));
    ws.on("message", (msg) => on_message.call(this, ws, msg));
  });

  // ping clients occasionally and terminate dead connections
  setup_ws_ping_interval.call(this);
}
