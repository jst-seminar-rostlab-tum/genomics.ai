import { BACKEND_ADDRESS, MULTIPART_UPLOAD_STATUS as Status, UPLOAD_CHUNK_SIZE } from '../utils/common/constants';
import removeItemFromArray from '../utils/common/utils';

export function getSubmissionProgressPercentage(progress) {
  if (progress.chunks === 0) {
    return 0;
  }
  return (progress.uploaded / progress.chunks) * 100;
}

function expectStatus(response, requestName, code) {
  if (response.status !== code) {
    throw new Error(`Invalid status for request "${requestName}" (${response.status})`);
  }
  return response;
}

function getAuthHeader() {
  return {
    auth: localStorage.getItem('jwt'),
  };
}

function getAuthAndJsonHeader() {
  return {
    auth: localStorage.getItem('jwt'),
    'content-type': 'application/json',
  };
}

async function uploadChunks(chunkCount, remaining, selectedFile, uploadId,
  submissionProgress, setSubmissionProgress, promiseArray) {
  for (let index = 1;
    index < chunkCount + 1 && !localStorage.getItem('cancelUpload');
    index += 1) {
    if (!remaining.includes(index)) {
      console.log(`Skipping chunk ${index}`);
      // eslint-disable-next-line no-continue
      continue;
    }
    const start = (index - 1) * UPLOAD_CHUNK_SIZE;
    const end = (index) * UPLOAD_CHUNK_SIZE;
    const blob = (index < chunkCount) ? selectedFile.slice(start, end)
      : selectedFile.slice(start);

    const currentPromise = fetch(`${BACKEND_ADDRESS}/file_upload/get_upload_url?${new URLSearchParams({
      partNumber: index,
      uploadId,
    })}`, { method: 'GET', headers: getAuthHeader() })
      .then((response) => expectStatus(response, 'get_upload_url', 200))
      .then((response) => response.json())
      .then(async ({ presignedUrl }) => fetch(presignedUrl, {
        body: blob,
        method: 'PUT',
        headers: { 'Content-Type': selectedFile.type },
      })
        .then((chunkResponse) => {
          if (chunkResponse.status !== 200) {
            throw new Error(`Invalid status code received for chunk upload (${chunkResponse.status})`);
          }
          const etag = chunkResponse.headers.get('ETag');
          const newRemaining = removeItemFromArray(submissionProgress.remaining, index);
          setSubmissionProgress((prevState) => (
            { ...prevState, uploaded: prevState.uploaded + 1, newRemaining }
          ));
          return { etag, index };
        }));
    promiseArray.push(currentPromise);
    try {
      // sequentializing upload
      // eslint-disable-next-line no-await-in-loop
      await currentPromise;
    } catch (ignored) { /* ignored */ }
  }
}

function finishUpload(chunkCount, promiseArray, submissionProgress, setSubmissionProgress,
  selectedFile, uploadId) {
  Promise.all(promiseArray.map((promise) => promise.catch((e) => e)))
    .then(async (promises) => {
      if (localStorage.getItem('cancelUpload')) {
        localStorage.removeItem('cancelUpload');
        return;
      }
      if (submissionProgress.remaining.length > 0) {
        setSubmissionProgress((prevState) => ({
          ...prevState,
          status: Status.ERROR_PROGRESS,
        }));
        return;
      }
      const uploadPartsArray = submissionProgress.uploadedParts;
      promises.forEach(({ etag, index }) => {
        uploadPartsArray.push({
          ETag: etag,
          PartNumber: index,
        });
      });
      setSubmissionProgress((prevState) => ({
        ...prevState, uploadedParts: uploadPartsArray, status: Status.UPLOAD_FINISHING,
      }));

      fetch(`${BACKEND_ADDRESS}/file_upload/complete_upload`, {
        method: 'POST',
        headers: getAuthAndJsonHeader(),
        body: JSON.stringify({
          parts: uploadPartsArray.sort((a, b) => (a.PartNumber - b.PartNumber)),
          uploadId,
        }),
      }).then((response) => expectStatus(response, 'complete_upload', 200))
        .then(() => setSubmissionProgress((prevState) => ({ ...prevState, status: 'complete' })))
        .catch(() => {
          setSubmissionProgress((prevState) => ({
            ...prevState,
            status: 'error_finish',
          }));
        });
    }).catch((err) => console.log(err));
}

export async function uploadMultipartFile(uploadId, selectedFile,
  submissionProgress, setSubmissionProgress) {
  const chunkCount = Math.floor(selectedFile.size / UPLOAD_CHUNK_SIZE) + 1;
  const promiseArray = [];

  setSubmissionProgress((prevState) => ({ ...prevState, chunks: chunkCount }));

  if (submissionProgress.status !== Status.ERROR_FINISH) {
    const { remaining } = submissionProgress;
    if (remaining.length === 0) {
      for (let i = 1; i < chunkCount + 1; i += 1) {
        remaining.push(i);
      }
      setSubmissionProgress((prevState) => ({ ...prevState, remaining }));
    }
    setSubmissionProgress((prevState) => ({ ...prevState, status: Status.UPLOAD_PROGRESS }));
    await uploadChunks(chunkCount, remaining, selectedFile,
      uploadId, submissionProgress, setSubmissionProgress, promiseArray);
  }
  finishUpload(chunkCount, promiseArray, submissionProgress,
    setSubmissionProgress, selectedFile, uploadId);
}

const startUpload = (selectedFile,
  submissionProgress,
  setSubmissionProgress,
  projectData = null) => fetch(
  `${BACKEND_ADDRESS}/file_upload/start_upload`,
  {
    method: 'POST',
    headers: getAuthAndJsonHeader(),
    body: JSON.stringify(projectData || { fileName: selectedFile.name }),
  },
)
  .then((response) => expectStatus(response, 'start_upload', 200))
  .then((response) => response.json())
  .then(({ uploadId }) => {
    setSubmissionProgress((prevState) => ({ ...prevState, uploadId }));
    return uploadMultipartFile(uploadId,
      selectedFile, submissionProgress, setSubmissionProgress);
  })
  .catch(() => {
    setSubmissionProgress((prevState) => ({
      ...prevState,
      status: Status.ERROR_START,
    }));
  });

export async function startOrContinueUpload(
  selectedFile,
  submissionProgress,
  setSubmissionProgress,
  projectData = null,
) {
  localStorage.removeItem('cancelUpload');
  const { status, uploadId } = submissionProgress;
  if (status === Status.ERROR_PROGRESS || status === Status.ERROR_FINISH) {
    return uploadMultipartFile(uploadId, selectedFile, submissionProgress, setSubmissionProgress);
  }
  return startUpload(selectedFile, submissionProgress, setSubmissionProgress, projectData);
}
