import ws from "ws";

declare module "ws" {
  export interface WebSocket extends ws {
    // Has this websocket timed out?
    // Every x seconds: Is set to true and "ping" is sent
    // When "pong" is received: Is set to false.
    timeout: boolean;
  }
}
