export enum WsMessageType {
    // from client - request to subscribe to a data channel
    WS_SUBSCRIBE = "subscribe",

    // from server - confirms subscription
    WS_SUBSCRIBED = "subscribed",

    // from client - request to unsubscribe from a data channel
    WS_UNSUBSCRIBE = "unsubscribe",

    // from server - confirms unsubscription
    WS_UNSUBSCRIBED = "unsubscribed",

    // from server - new data for subscribed data channel
    WS_DATA = "data",

    // from any - an error has occurred
    WS_ERROR = "error",
}

export type WsChannelData = {
    channel : string,
    data : any,
}

export type WsMessage = {
    type : WsMessageType,
    data : any
};

export function parseWsMessage(text : string){
    let obj = JSON.parse(text);
    if(!obj.hasOwnProperty("type"))
        throw new Error( "Error parsing WsMessage: missing 'type' property");
    if(!Object.values(WsMessageType).includes(obj.type))
        throw new Error( "Error parsing WsMessage: 'type' property invalid");
    if(!obj.hasOwnProperty("data"))
        throw new Error( "Error parsing WsMessage: missing 'data' property");

    let message : WsMessage = {
        type: obj.type,
        data: obj.data
    };
    return message;
}