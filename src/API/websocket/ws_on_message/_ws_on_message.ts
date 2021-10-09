import Websocket_Host from "../_websocket_host";
import WebSocket from "ws";
import {parseWsMessage, WsMessage, WsMessageType} from "../types/ws_message";

// event handler for receiving a new message from a connected websocket
export function on_message(this:Websocket_Host, ws:WebSocket, message: WebSocket.Data){
    // ******* Parse message *******
    let ws_msg : WsMessage;
    try{
        ws_msg = parseWsMessage(message.toString());
    }catch(e){
        return ws.send(JSON.stringify({
            type: "error",
            message: "unable to parse json-data"
        }));
    }

    // ******* Process message *******
    try{
        switch (ws_msg.type){
            default:
                this.ws_send(ws, WsMessageType.WS_ERROR, {message: "invalid message type '" + ws_msg.type + "'"});
        }
    }catch(e){
        this.ws_send(ws, WsMessageType.WS_ERROR, {message: "failed to process message"});
    }
}