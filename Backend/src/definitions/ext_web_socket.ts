import WebSocket from "ws";

export interface ExtWebSocket extends WebSocket {
    timeout?: boolean; // declare optional property "timeout"
}