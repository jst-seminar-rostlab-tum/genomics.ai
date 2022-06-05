import { HTTP_Server } from "../../http_server/http_server";

import express from "express";
import * as core from "express-serve-static-core";
import { express_middleware, express_routes_middleware } from "./express_middleware";
import { express_routes } from "./express_routes";

export default class REST_Host {
  // ---- fields ----
  // instance of the http-server hosting the websocket connections
  protected http_server: HTTP_Server;
  // Express server-instance
  protected expressApp: core.Express;
  // ------------

  constructor(http_server: HTTP_Server) {
    this.http_server = http_server;
    this.expressApp = express();
  }

  // performs all necessary setup
  public async init() {
    this.expressApp.use(express_middleware.call(this));
    this.expressApp.use("/v1", express_routes_middleware.call(this), express_routes() );

    // forward http-requests to express.js server
    this.http_server.server?.on("request", this.expressApp);

    console.log(`REST-Server operational.`);
    return this;
  }
}
