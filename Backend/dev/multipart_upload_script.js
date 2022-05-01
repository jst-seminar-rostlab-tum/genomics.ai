//********************************************//
// Example of how to use the multipart routes* //
//********************************************//

const https = require("https");
const http = require("http");
const fs = require("fs");
const querystring = require("querystring");
const argsJson = require("./scriptargs.json");

const HOST = "http://localhost:8050";

let email = argsJson.email;
let pwd = argsJson.password;
let filepath = argsJson.filepath;
let projectName = argsJson.projectName;
let atlasId = argsJson.atlasId;
let modelId = argsJson.modelId;
let chunkSize = Number.parseInt(argsJson.chunkSize);

console.log("Loading file");
let fileData = fs.readFileSync(filepath);
console.log("Loaded file");
let chunkNum = Math.ceil(fileData.byteLength / chunkSize);
console.log(chunkNum + " chunks in file");

do_upload();

function request(method, path, jwt, contentType, body, query) {
  return new Promise((resolve, reject) => {
    if (query) {
      path += "?" + new URLSearchParams(query);
    }
    let headers = {};
    if (contentType) headers["Content-Type"] = contentType;
    let rawBody;
    let contentLength;
    if (body) {
      if (body instanceof Buffer || body instanceof Uint8Array || body instanceof String) {
        rawBody = body;
      } else {
        rawBody = JSON.stringify(body);
      }
      contentLength = rawBody.length;
      headers["Content-Length"] = contentLength;
    }
    if (jwt) {
      headers.auth = jwt;
    }
    let module = path.startsWith("https") ? https : http;
    let req = module.request(
      path,
      {
        method,
        headers: headers,
      },
      function (res) {
        let data = "";
        res.on("data", (chunk) => {
          data += chunk;
        });
        res.on("end", () => {
          resolve({ body: data, headers: res.headers });
        });
      }
    );
    if (body) {
      req.write(rawBody);
    }
    req.end();
  });
}

async function do_upload() {
  let jwt = await auth(email, pwd);

  let project = await start_upload(jwt, projectName, filepath, atlasId, modelId);
  let uploadId = project.uploadId;

  let parts = [];

  let randomOrder = [];
  for (let i = 1; i <= chunkNum; i++) {
    randomOrder.push(i);
  }
  shuffle(randomOrder);
  console.log(randomOrder);

  for (let i = 0; i < chunkNum; i++) {
    let partNum = randomOrder[i];
    let from = (partNum - 1) * chunkSize;
    let to = Math.min(partNum * chunkSize, fileData.byteLength);
    let chunk = fileData.slice(from, to);
    let presignedUrl = await get_upload_url(jwt, uploadId, partNum);
    let etag = await put_to_presigned_url(presignedUrl, chunk);
    parts.push({ PartNumber: partNum, ETag: etag });
  }
  parts.sort((a, b) => a.PartNumber - b.PartNumber);
  await complete_upload(jwt, uploadId, parts);
  while (true) {
    let projectUpdated = await get_project_info(jwt, project._id);
    if (projectUpdated.status === "PROCESSING_PENDING") {
      console.log("Not done yet");
      await wait(4000);
    } else if (projectUpdated.status === "DONE") {
      console.log("PROCESSING IS DONE => get result");
      await get_from_presigned_url(projectUpdated.location);
      break;
    } else {
      console.log("Wrong state");
      break;
    }
  }
}

async function wait(ms) {
  await new Promise((resolve, rej) => setTimeout(resolve, ms));
}

async function auth(email, pwd) {
  console.log("===================================");
  console.log("Sending auth to " + HOST + "/auth");
  let { body } = await request("POST", HOST + "/auth", null, "application/json", {
    email: email,
    password: pwd,
  });
  try {
    let { jwt } = JSON.parse(body);
    console.log("Received jwt token");
    return jwt;
  } catch (e) {
    console.error("Error while parsing response");
    throw e;
  }
}

async function start_upload(jwt, projectName, filepath, atlasId, modelId) {
  console.log("===================================");
  console.log("Starting upload");
  let { body } = await request(
    "POST",
    HOST + "/file_upload/start_upload",
    jwt,
    "application/json",
    {
      projectName,
      fileName: filepath,
      atlasId: atlasId,
      modelId: modelId,
    }
  );
  console.log("Received response:");
  try {
    console.log(JSON.parse(body));
  } catch {
    console.log(body);
  }
  try {
    let project = JSON.parse(body);
    console.log("UploadId: " + project.uploadId);
    return project;
  } catch (e) {
    console.error("Error while parsing response");
    throw e;
  }
}

async function get_upload_url(jwt, uploadId, partNumber) {
  console.log("===================================");
  console.log("Getting url for part " + partNumber);
  let { body } = await request(
    "GET",
    HOST + "/file_upload/get_upload_url",
    jwt,
    undefined,
    undefined,
    {
      uploadId,
      partNumber,
    }
  );
  console.log("Received response:");
  try {
    console.log(JSON.parse(body));
  } catch {
    console.log(body);
  }
  try {
    let { presignedUrl } = JSON.parse(body);
    console.log("URL: " + presignedUrl);
    return presignedUrl;
  } catch (e) {
    console.error("Error while parsing response");
    throw e;
  }
}

async function complete_upload(jwt, uploadId, parts) {
  console.log("===================================");
  console.log("Completing upload: Parts: ");
  console.log(parts);
  let { body } = await request(
    "POST",
    HOST + "/file_upload/complete_upload",
    jwt,
    "application/json",
    {
      uploadId,
      parts,
    }
  );
  console.log("Received response:");
  console.log(body);
}

async function put_to_presigned_url(url, buffer) {
  console.log("===================================");
  console.log("Uploading to presigned url");
  let { body, headers } = await request("PUT", url, undefined, "application/octet-stream", buffer);
  console.log("Received response");
  console.log("Headers:");
  console.log(headers);
  console.log("Body:");
  console.log(body);
  return headers.etag;
}

async function get_from_presigned_url(url) {
  console.log("=================================");
  console.log("Downloading Result file...");
  let { body } = await request("GET", url);
  if (body.length < 1000) {
    console.log(body);
  } else {
    console.log("<Large file>");
  }
}

async function get_project_info(jwt, id) {
  console.log("=================================");
  console.log("Getting Project state...");
  let { body } = await request("GET", HOST + "/project/" + id, jwt);

  try {
    console.log(JSON.parse(body));
  } catch {
    console.log(body);
  }

  try {
    let project = JSON.parse(body);
    console.log("STATUS: " + project.status);
    return project;
  } catch (e) {
    console.error("Error while parsing response");
    throw e;
  }
}

function shuffle(a) {
  //Fisher-Yates
  for (let i = a.length - 1; i >= 1; i--) {
    let j = Math.floor(Math.random() * (i + 1));
    let tmp = a[j];
    a[j] = a[i];
    a[i] = tmp;
  }
}
