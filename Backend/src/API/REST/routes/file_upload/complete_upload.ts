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

export default function upload_complete_upload_route() {
  let router = express.Router();
  router.post("/file_upload/complete_upload", check_auth(), async (req: ExtRequest, res) => {
    let { parts, uploadId } = req.body;
    if (!process.env.S3_BUCKET_NAME) return res.status(500).send("Server was not set up correctly");

    try {
      let project = await ProjectService.getProjectByUploadId(String(uploadId), req.user_id);

      if (!project) return res.status(400).send("Project could not be found");

      //Complete multipart upload
      let params: CompleteMultipartUploadRequest = {
        Bucket: process.env.S3_BUCKET_NAME,
        Key: `projects/${project._id}/query.h5ad`,
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
        let [model, atlas] = await Promise.all([
          ModelService.getModelById(project.modelId),
          AtlasService.getAtlasById(project.atlasId),
        ]);
        if (!model) {
          return res.status(500).send("Could not find model");
        }
        let queryInfo = {
          model: model.name,
          query_data: `projects/${project.id}/query.h5ad`,
          output_path: `results/${project.id}/query.tsv`,
          model_path: `results/${project.id}/model.pt`,
          reference_data: `atlas/${project.atlasId}/data.h5ad`,
          //ref_path: `models/${project.modelId}/model.pt`,
          async: false,
        };
        if (process.env.CLOUD_RUN_URL) {
          const url = `${process.env.CLOUD_RUN_URL}/query`;
          const auth = new GoogleAuth();
          const client = await auth.getIdTokenClient(url);
          //Send response before processing
          res.status(200).send("Processing started");
          //Processing is synchronous, response is sent by ML only after the result is produced, might take some time
          await client.request({ url, method: "POST", body: JSON.stringify(queryInfo) });
        } else if (process.env.NODE_ENV!="production") {
          console.log(
            "CLOUD_RUN_URL not defined, falling back to dummy result with 10s processing time"
          );
          res.status(200).send("Processing started");
          const params2: PutObjectRequest = {
            Bucket: process.env.S3_BUCKET_NAME!,
            Key: `results/${project!.id}/query.tsv`,
            Body: "This is a dummy result",
          };
          await s3.upload(params2).promise();
          await new Promise((resolve, reject) => {
            setTimeout(resolve, 10000);
          });
        } else {
          const updateStatus: UpdateProjectDTO = { status: ProjectStatus.PROCESSING_FAILED };
          await ProjectService.updateProjectByUploadId(params.UploadId, updateStatus);
          return res.status(500).send("Processing failed!");
        }

        //Processing finished, http response has already be sent before processing, update database entry now
        let params2: any = {
          Bucket: process.env.S3_BUCKET_NAME!,
          Key: `results/${project!.id}/query.tsv`,
          Expires: 60 * 60 * 24 * 7 - 1, // one week minus one second
        };
        let presignedUrl = await s3.getSignedUrlPromise("getObject", params2);
        const updateLocationAndStatus: UpdateProjectDTO = {
          location: presignedUrl,
          status: ProjectStatus.DONE,
        };
        await ProjectService.updateProjectByUploadId(params.UploadId, updateLocationAndStatus);
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
  });
  return router;
}
