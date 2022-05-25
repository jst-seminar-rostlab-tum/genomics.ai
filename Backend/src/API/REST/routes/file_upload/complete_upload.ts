import express from "express";
import {
  CompleteMultipartUploadOutput,
  CompleteMultipartUploadRequest,
  PutObjectRequest,
} from "aws-sdk/clients/s3";
import { AWSError, S3 } from "aws-sdk";

import check_auth from "../../middleware/check_auth";
import { ExtRequest } from "../../../../definitions/ext_request";
import ProjectService from "../../../../database/services/project.service";
import { UpdateProjectDTO } from "../../../../database/dtos/project.dto";
import s3 from "../../../../util/s3";
import { GoogleAuth } from "google-auth-library";
import { ProjectStatus } from "../../../../database/models/project";
import AtlasService from "../../../../database/services/atlas.service";
import ModelService from "../../../../database/services/model.service";

import { validationMdw } from "../../middleware/validation";
import ProjectUpdateTokenService from "../../../../database/services/project_update_token.service";
import { query_path, result_model_path, result_path } from "./bucket_filepaths";

export default function upload_complete_upload_route() {
  let router = express.Router();
  router.post(
    "/file_upload/complete_upload",
    validationMdw,
    check_auth(),
    async (req: ExtRequest, res) => {
      let { parts, uploadId } = req.body;
      if (!process.env.S3_BUCKET_NAME)
        return res.status(500).send("Server was not set up correctly");

      try {
        let project = await ProjectService.getProjectByUploadId(String(uploadId), req.user_id);

        if (!project) return res.status(400).send("Project could not be found");

        //Complete multipart upload
        let params: CompleteMultipartUploadRequest = {
          Bucket: process.env.S3_BUCKET_NAME,
          Key: query_path(project._id),
          MultipartUpload: { Parts: parts },
          UploadId: String(uploadId),
        };
        let data;
        try {
          data = await s3.completeMultipartUpload(params).promise();
        } catch (err: any) {
          console.error(err, err.stack || "Error when completing multipart upload");
          return res.status(500).send(err);
        }
        if (!data || !data.Key || !data.Bucket || !data.Location)
          return res.status(500).send("Error getting Multipart-Upload object data");

        //Query file size and save in project
        try {
          let request: S3.HeadObjectRequest = { Key: data.Key, Bucket: data.Bucket };
          let result = await s3.headObject(request).promise();
          const updateFileAndStatus: UpdateProjectDTO = {
            fileSize: result.ContentLength,
            status: ProjectStatus.PROCESSING_PENDING,
          };
          await ProjectService.updateProjectByUploadId(params.UploadId, updateFileAndStatus);
          if (process.env.CLOUD_RUN_URL) {
            let [model, atlas] = await Promise.all([
              ModelService.getModelById(project.modelId),
              AtlasService.getAtlasById(project.atlasId),
            ]);
            if (!model || !atlas) {
              await ProjectService.updateProjectById(params.UploadId, {
                status: ProjectStatus.PROCESSING_FAILED,
              });
              return res.status(500).send(`Could not find ${!model ? "model" : "atlas"}`);
            }

            //Create a token, which can be used later to update the projects status
            let { token: updateToken } = await ProjectUpdateTokenService.addToken({
              _projectId: project._id,
            });
            let queryInfo = {
              model: model.name,
              query_data: query_path(project.id),
              output_path: result_path(project.id),
              model_path: result_model_path(project.id),
              reference_data: `atlas/${project.atlasId}/data.h5ad`,
              //ref_path: `models/${project.modelId}/model.pt`,
              async: false,
              webhook: `${process.env.API_URL}/projects/results/${updateToken}`,
            };
            console.log("sending: ");
            console.log(queryInfo);
            const url = `${process.env.CLOUD_RUN_URL}/query`;
            const auth = new GoogleAuth();
            const client = await auth.getIdTokenClient(url);
            //Send response before processing
            res.status(200).send("Processing started");
            //Processing is synchronous, response is sent by ML only after the result is produced, might take some time
            let result;
            try {
              result = await client.request({
                url,
                method: "POST",
                body: JSON.stringify(queryInfo),
              });
            } catch (e: any) {
              console.log("Processing failed:");
              console.log(e.message || e);
              result = null;
            }
            if (!result || result.status != 200) {
              await ProjectService.updateProjectByUploadId(params.UploadId, {
                status: ProjectStatus.PROCESSING_FAILED,
              });
              return;
            }
          } else if (process.env.NODE_ENV != "production") {
            console.log(
              "CLOUD_RUN_URL not defined, falling back to dummy result with 10s processing time"
            );
            res.status(200).send("Processing started");

            try {
              let fs = await import("fs");
              let path = await import("path");
              let content: Buffer = await new Promise((resolve, reject) => {
                fs.readFile(
                  path.join(__dirname, "../../../../../dev/test_file1.csv"),
                  function (err, data) {
                    if (err) reject(err);
                    else resolve(data);
                  }
                );
              });
              const params2: PutObjectRequest = {
                Bucket: process.env.S3_BUCKET_NAME!,
                Key: result_path(project.id),
                Body: content,
              };
              await s3.upload(params2).promise();
              await new Promise((resolve, reject) => {
                setTimeout(resolve, 10000);
              });
            } catch (e) {
              await ProjectService.updateProjectByUploadId(params.UploadId, {
                status: ProjectStatus.PROCESSING_FAILED,
              });
              return;
            }
          } else {
            const updateStatus: UpdateProjectDTO = { status: ProjectStatus.PROCESSING_FAILED };
            await ProjectService.updateProjectByUploadId(params.UploadId, updateStatus);
            return res.status(500).send("Processing failed!");
          }
        } catch (err) {
          console.log(err);
          try {
            res.status(500).send(`Error persisting Multipart-Upload object data: ${err}`);
          } catch {}
        }
      } catch (err) {
        console.log(err);
        res.status(500).send(err);
      }
    }
  );
  return router;
}
