import Ajv, { JSONSchemaType } from "ajv";
export const ajv = new Ajv();

interface Auth {
  email: string;
  password: string;
}
const authSchema: JSONSchemaType<Auth> = {
  type: "object",
  properties: {
    email: { type: "string" },
    password: { type: "string" },
  },
  required: ["email", "password"],
  additionalProperties: false,
};

ajv.addSchema(authSchema, "auth").getSchema("auth");
