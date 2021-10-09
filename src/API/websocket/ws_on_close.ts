import Websocket_Host from "./_websocket_host";
import WebSocket from "ws";

// event handler for when a connected websocket closes the connection
export function on_close(this:Websocket_Host, ws:WebSocket){
    // unsubscribe socket on close
    for (let channel in this.subscriptions) {
        // remove socket from "subscribers" array where it exists
        if(this.subscriptions[channel].subscribers.includes(ws))
            this.subscriptions[channel].subscribers.splice(
                this.subscriptions[channel].subscribers.indexOf(ws), 1);
    }
}