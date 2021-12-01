// export const BACKEND_ADDRESS = 'https://custom-helix-329116.ey.r.appspot.com';
export const BACKEND_ADDRESS = 'http://localhost:8050';

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
};
export const statusIsUpload = (status) => status.startsWith('upload_');
export const statusIsError = (status) => status.startsWith('error_');
