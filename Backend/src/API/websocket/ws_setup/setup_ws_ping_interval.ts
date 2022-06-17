import Websocket_Host from "../_websocket_host";
import { ExtWebSocket } from "../../../definitions/ext_web_socket";

// setup the periodic pinging of all connected websockets
export function setup_ws_ping_interval(this: Websocket_Host) {
  const self = this;

  this.ws_ping_interval = setInterval(() => {
    self.wss?.clients.forEach((ws: ExtWebSocket) => {
      // terminate timed out sockets (= sent ping last interval and received no pong yet)
      if (ws.timeout) return ws.terminate();

      // Send Ping and set ws.timeout to "true" (ws.timeout is reset to "false" when pong is received)
      ws.timeout = true;
      ws.ping();
    });
  }, 30 * 1000); // every 30 secs

  this.wss?.on("close", () => {
    if (self.ws_ping_interval) clearInterval(self.ws_ping_interval);
  });
}
