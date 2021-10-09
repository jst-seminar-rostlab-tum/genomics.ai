import Websocket_Host from "./_websocket_host";
import {WsMessage, WsMessageType} from "./types/ws_message";
import WebSocket from "ws";

// send a data-object of a specific message type to a connected client
export function ws_send(this:Websocket_Host, ws:WebSocket, type:WsMessageType, data:any){
    let msg : WsMessage = {
        type: type,
        data: data
    };

    return ws.send(JSON.stringify(msg));
}