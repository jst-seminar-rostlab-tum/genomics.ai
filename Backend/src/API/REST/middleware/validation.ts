import Ajv from "ajv";
import { Request, Response, NextFunction } from "express";
import * as swaggerDocument from "../../../swagger.json";

const ajv = new Ajv();

for (const [path, pathObj] of Object.entries(swaggerDocument.paths)) {
  for (const [method, methodObj] of Object.entries(pathObj)) {
    const schemaName = constructSchemaName(path, method);
    const schema = methodObj.requestBody.content["application/json"].schema;
    ajv.addSchema(schema, schemaName).getSchema(schemaName);
  }
}

export function validationMdw(req: Request, res: Response, next: NextFunction) {
  const schemaName = constructSchemaName(req.path, req.method);

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
  return `${path}_${method}`.toLowerCase();
}
