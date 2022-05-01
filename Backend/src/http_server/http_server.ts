import http from "http";
import https from "https";
import fs from "fs/promises";

export class HTTP_Server {
  server: http.Server | undefined = undefined;

  // setup the http(s)-server which hosts the websocket-server and REST-APIs
  async setup() {
    if (
      ["SSL_CERT", "SSL_KEY", "SSL_CA"].map((x) => process.env.hasOwnProperty(x)).includes(false)
    ) {
      if (process.env.NODE_ENV="production") throw new Error("SSL-env-variables missing!");
      else {
        console.log("Running unencrypted HTTP-server! (dev mode)");
        this.server = http.createServer();
      }
    } else {
      const privateKey = await fs.readFile(process.env.SSL_KEY + "", "utf8");
      const certificate = await fs.readFile(process.env.SSL_CERT + "", "utf8");
      const ca = await fs.readFile(process.env.SSL_CA + "", "utf8");

      this.server = https.createServer({
        key: privateKey,
        cert: certificate,
        ca: ca,
      });
    }

    this.server.listen({ port: parseInt("" + process.env.PORT) }, () =>
      console.log("Server listening on port " + process.env.PORT)
    );
  }
}
