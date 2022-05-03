import * as React from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import {
  Box, IconButton, LinearProgress, Stack,
} from '@mui/material';
import CircleIcon from '@mui/icons-material/Circle';
import { useHistory } from 'react-router-dom';
import { getSubmissionProgressPercentage } from 'shared/services/UploadLogic';
import {
  MULTIPART_UPLOAD_STATUS, MULTIPART_UPLOAD_STATUS as Status, statusIsError, statusIsUpload, PROJECT_STATUS,
} from 'shared/utils/common/constants';
import Clear from '@mui/icons-material/Clear';
import ReplayIcon from '@mui/icons-material/Replay';
import ProgressBar from 'components/ProgressBar';

function ProcessingStatus() {
  return (
    <>
      <Box sx={{ pr: 2, flexGrow: 1 }}>
        <LinearProgress />
      </Box>
      <Typography variant="caption" noWrap>Processing by scArches...</Typography>
    </>
  );
}

export default function ProjectBarCard({
  projectId, name, status, submissionProgress, setSubmissionProgress,
}) {
  const history = useHistory();

  const color = status === PROJECT_STATUS.DONE
    ? 'lightGreen'
    : status === PROJECT_STATUS.ABORTED
    || (!submissionProgress && status === PROJECT_STATUS.UPLOAD_PENDING)
    || status === PROJECT_STATUS.PROCESSING_FAILED
      ? 'red'
      : 'orange';

  return (
    <Card sx={{
      marginTop: '5', marginBottom: '0.5em', borderStyle: 'solid', borderColor: '#C8C8C8', borderWidth: '0.1px',
    }}
    >
      <Stack direction="row" sx={{ justifyContent: 'space-between' }}>
        <Stack
          direction="row"
          spacing={4}
          sx={{ alignItems: 'center', flexGrow: 1 }}
        >

          <CircleIcon sx={{
            fontSize: 30, marginLeft: '3%', color,
          }}
          />
          <Typography noWrap sx={{ width: '30%' }}>
            {name}
          </Typography>
          {submissionProgress ? (
            <Box
              sx={{
                flexGrow: 1, display: 'flex', alignItems: 'center',
              }}
            >
              {statusIsUpload(submissionProgress.status)
                  && (
                  <>
                    <Box sx={{ pr: 2, flexGrow: 1 }}>
                      <LinearProgress variant="determinate" value={getSubmissionProgressPercentage(submissionProgress)} />
                    </Box>
                    <Typography variant="caption">Uploading...</Typography>
                    <IconButton
                      onClick={() => {
                        setSubmissionProgress((prevState) => (
                          { ...prevState, status: Status.CANCELING }));
                        localStorage.setItem('cancelUpload', '1'); // worst design ever
                      }}
                    >
                      <Clear color="error" />
                    </IconButton>
                  </>
                  )}
              {statusIsError(submissionProgress.status)
              && (
              <>
                <Typography>{submissionProgress.status}</Typography>
                <IconButton onClick={() => {}}>
                  <ReplayIcon />
                </IconButton>
              </>
              )}
              {submissionProgress.status === MULTIPART_UPLOAD_STATUS.COMPLETE
                   && <ProcessingStatus />}
            </Box>
          ) : null}
          {!submissionProgress
            ? (
              <Box
                sx={{ flexGrow: 1, display: 'flex', alignItems: 'center' }}
              >
                {status === PROJECT_STATUS.UPLOAD_PENDING
                   && (
                   <Typography variant="caption">
                     Upload failed or canceled
                   </Typography>
                   )}
                {status === PROJECT_STATUS.PROCESSING_PENDING
                   && <ProcessingStatus />}
                {(status === PROJECT_STATUS.ABORTED || status === PROJECT_STATUS.PROCESSING_FAILED)
                   && <Typography variant="caption">Processing failed</Typography>}
              </Box>
            )
            : null}
        </Stack>
        <Box sx={{
          p: 0.1, bgcolor: 'background.paper', borderRadius: 3, width: 'flex', mr: 3, display: 'flex', m: 2, alignItems: 'center',
        }}
        >
          <Button
            variant="outlined"
            // size="small"
            sx={{
              borderRadius: 100,
              mr: 2,
            }}
            style={{ textTransform: 'none' }}
          >
            Add To Team
          </Button>
          <Button
            variant="contained"
            color="secondary"
            sx={{
              borderRadius: 100,
            }}
            style={{ textTransform: 'none' }}
            disabled={status !== 'DONE'}
            onClick={() => history.push(`./genemapper/result/${projectId}`)}
          >
            See Results
          </Button>
        </Box>
      </Stack>
    </Card>
  );
}
