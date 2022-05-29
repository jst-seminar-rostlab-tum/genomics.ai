import {
  red,
  green,
  yellow,
  grey,
  blue,
} from '@mui/material/colors';

export const BACKEND_ADDRESS = 'https://devversion-dot-custom-helix-329116.ey.r.appspot.com/v1';

export const JOB_QUEUE_UPDATE_INTERVAL = 5000;
export const PROJECTS_UPDATE_INTERVAL = 5000;

export const UPLOAD_CHUNK_SIZE = 25000000; // 50MB
export const MULTIPART_UPLOAD_STATUS = {
  IDLE: 'idle',
  SELECTED: 'selected',
  UPLOAD_STARTING: 'upload_starting',
  ERROR_START: 'error_start',
  UPLOAD_PROGRESS: 'upload_progress',
  ERROR_PROGRESS: 'error_progress',
  UPLOAD_FINISHING: 'upload_finishing',
  ERROR_FINISH: 'error_finish',
  COMPLETE: 'complete',
  CANCELING: 'canceling',
};
export const statusIsUpload = (status) => status.startsWith('upload_');
export const statusIsError = (status) => status.startsWith('error_');

export const jobStatusColors = {
  UPLOAD_PENDING: yellow[600],
  PROCESSING_PENDING: blue[300],
  ABORTED: red[300],
  DONE: green[300],
  unknown: grey[500],
};

export const jobStatusTitles = {
  UPLOAD_PENDING: 'UPLOADING',
  UPLOAD_COMPLETE: 'UPLOADED',
  PROCESSING_PENDING: 'PROCESSING',
  ABORTED: 'CANCELLED',
  DONE: 'COMPLETED',
  unknown: 'UNKNOWN',
};

export const PROJECT_STATUS = {
  UPLOAD_PENDING: 'UPLOAD_PENDING',
  UPLOAD_COMPLETE: 'UPLOAD_COMPLETE',
  PROCESSING_PENDING: 'PROCESSING_PENDING',
  PROCESSING_FAILED: 'PROCESSING_FAILED',
  ABORTED: 'ABORTED',
  DONE: 'DONE',
};
