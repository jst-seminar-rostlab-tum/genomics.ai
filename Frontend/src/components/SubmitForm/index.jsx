import React, { useState } from 'react';
import LoadingButton from '@mui/lab/LoadingButton';
import UploadIcon from '@mui/icons-material/Upload';
import ReplayIcon from '@mui/icons-material/Replay';
import { Typography } from '@mui/material';
import { Clear } from '@mui/icons-material';
import styles from './submitform.module.css';
import UploadField from '../UploadField';
import ProgressBar from '../ProgressBar';
import { getSubmissionProgressPercentage, startOrContinueUpload } from 'shared/services/UploadLogic';
import { MULTIPART_UPLOAD_STATUS as Status, statusIsError, statusIsUpload } from 'shared/utils/common/constants';

function SubmitForm() {
  const [selectedFile, setSelectedFile] = useState();
  const [submissionProgress, setSubmissionProgress] = useState({
    status: Status.IDLE,
    uploadId: '',
    chunks: 0,
    uploaded: 0,
    remaining: [],
    uploadedParts: [],
  });

  const changeHandler = (event) => {
    if (event.target.files[0]) {
      setSubmissionProgress({
        status: Status.SELECTED,
        uploadId: '',
        chunks: 0,
        uploaded: 0,
        remaining: [],
        uploadedParts: [],
      });
      setSelectedFile(event.target.files[0]);
    }
  };

  const {
    status, uploaded, chunks, remaining,
  } = submissionProgress;

  return (
    <div>
      <UploadField
        disabled={statusIsUpload(status) || status === Status.CANCELING}
        selectedFile={selectedFile}
        onChange={changeHandler}
      />
      <Typography>
        {
          {
            idle: 'Select a file to upload!',
            selected: 'Press Upload to start uploading!',
            upload_starting: 'Starting upload...',
            error_start: 'Couldn\'t start upload :/',
            upload_progress: `Uploading chunks... (${uploaded} of ${chunks})`,
            error_progress: `${remaining.length} chunks couldn't be uploaded!`,
            upload_finishing: 'Processing upload...',
            error_finish: 'Couldn\'t finish upload :/',
            canceling: 'Canceling upload...',
            complete: 'Upload complete!',
          }[status]
        }
      </Typography>
      <div className={styles.flexContainer}>
        {status !== 'idle' && status !== 'selected' ? <ProgressBar value={getSubmissionProgressPercentage(submissionProgress)} /> : null}
      </div>
      <div>
        <LoadingButton
          disabled={status === Status.IDLE || status === Status.COMPLETE}
          loading={statusIsUpload(status)}
          color={statusIsError(status) ? 'error' : 'primary'}
          startIcon={statusIsError(status) ? <ReplayIcon /> : <UploadIcon />}
          onClick={() => {
            setSubmissionProgress((prevState) => (
              { ...prevState, status: Status.UPLOAD_STARTING }));
            return startOrContinueUpload(selectedFile, submissionProgress, setSubmissionProgress);
          }}
        >
          {statusIsError(status) ? 'Retry' : 'Upload'}
        </LoadingButton>
        {((statusIsUpload(status) && status !== Status.UPLOAD_FINISHING)
          || status === Status.CANCELING)
          && (
            <LoadingButton
              disabled={status === Status.CANCELING}
              loading={status === Status.CANCELING}
              color="error"
              startIcon={<Clear />}
              onClick={() => {
                setSubmissionProgress((prevState) => (
                  { ...prevState, status: Status.CANCELING }));
                localStorage.setItem('cancelUpload', '1'); // worst design ever
                setSelectedFile(null);
              }}
            >
              Cancel
            </LoadingButton>
          )}
      </div>
    </div>
  );
}

export default SubmitForm;
