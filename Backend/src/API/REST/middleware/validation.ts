import Ajv from "ajv";
import { Request, Response, NextFunction } from "express";

import { loadSwaggerDocument } from "../../../swagger/load-swagger";

const ajv = new Ajv();

const swaggerDocument = loadSwaggerDocument();

for (const [path, pathObj] of Object.entries(swaggerDocument.paths)) {
  for (const [method, methodObj] of Object.entries(pathObj)) {
    const schemaName = constructSchemaName(path, method);

    const reqBody = (methodObj as any).requestBody;

    if (!reqBody) continue;
    let content = reqBody.content;
    if (!content && reqBody.$ref) {
      content = resolveRef(swaggerDocument, reqBody.$ref);
    }
    let schema = content?.["application/json"]?.schema;
    if (schema) ajv.addSchema(schema, schemaName).getSchema(schemaName);
  }
}

export function validationMdw(req: Request, res: Response, next: NextFunction) {
  const schemaName = constructSchemaName(req.route.path, req.method);

  const validate = ajv.getSchema(schemaName);

  if (!validate) console.warn("Schema not defined for", schemaName);

  if (validate && !validate(req.body)) {
    console.error(validate?.errors);
    return res.status(400).send({
      errors: validate?.errors,
    });
  }
  next();
}

function constructSchemaName(path: string, method: string) {
  path = path.replace("{", ":");
  path = path.replace("}", "");
  return `${path}_${method}`.toLowerCase();
}

function resolveRef(doc: any, ref: string): any {
  if (!ref.startsWith("#/")) {
    throw new Error("invalid ref");
  }
  let obj = doc;
  let components = ref.substr(2).split("/");
  for (const component of components) {
    obj = obj[component];
  }
  return obj;
}
