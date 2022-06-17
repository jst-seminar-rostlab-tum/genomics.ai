const fs = require("fs");
const path = require("path");

const SWAGGER_DIR = "../../swagger";
const SWAGGER_PATHS_DIR = `${SWAGGER_DIR}/paths`;

let swaggerMainDoc = null;

export function loadSwaggerDocument() {
  if (swaggerMainDoc) return swaggerMainDoc;

  let filePath: string | undefined;
  let pathContent: any;
  try {
    swaggerMainDoc = getSwaggerMainDoc();

    for (const p of fs.readdirSync(path.join(__dirname, SWAGGER_PATHS_DIR))) {
      filePath = path.join(__dirname, SWAGGER_PATHS_DIR, p);
      pathContent = JSON.parse(fs.readFileSync(filePath, "utf-8"));
      swaggerMainDoc.paths = { ...pathContent, ...swaggerMainDoc.paths };
    }
    return swaggerMainDoc;
  } catch (error) {
    console.error("trying to read the filePath...", filePath);
    throw error;
  }
}

function getSwaggerMainDoc() {
  const filePath = path.join(__dirname, `${SWAGGER_DIR}/swagger.json`);
  return JSON.parse(fs.readFileSync(filePath, "utf-8"));
}
