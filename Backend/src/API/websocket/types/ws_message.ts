export enum WsMessageType {
  // from client - request to authenticate client with JWT
  WS_AUTH = "authenticate",

  // from server - confirms authentication
  WS_AUTHED = "authenticated",

  // authentication was denied after WS_AUTH received
  WS_DENIED = "auth_denied",

  // from any - an error has occurred
  WS_ERROR = "error",
}

export type WsChannelData = {
  channel: string;
  data: any;
};

export type WsMessage = {
  type: WsMessageType;
  data: any;
};

export function parseWsMessage(text: string) {
  let obj = JSON.parse(text);
  if (!obj.hasOwnProperty("type"))
    throw new Error("Error parsing WsMessage: missing 'type' property");
  if (!Object.values(WsMessageType).includes(obj.type))
    throw new Error("Error parsing WsMessage: 'type' property invalid");
  if (!obj.hasOwnProperty("data"))
    throw new Error("Error parsing WsMessage: missing 'data' property");

  let message: WsMessage = {
    type: obj.type,
    data: obj.data,
  };
  return message;
}
