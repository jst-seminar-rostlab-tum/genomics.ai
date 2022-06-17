import Websocket_Host from "./_websocket_host";
import WebSocket from "ws";

// event handler for when a connected websocket closes the connection
export function on_close(this: Websocket_Host, ws: WebSocket) {}
