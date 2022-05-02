import * as React from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import {
  Box, CircularProgress, IconButton, LinearProgress, Stack,
} from '@mui/material';
import CircleIcon from '@mui/icons-material/Circle';
import { useHistory } from 'react-router-dom';
import { getSubmissionProgressPercentage } from 'shared/services/UploadLogic';
import CancelIcon from '@mui/icons-material/Cancel';
import {
  MULTIPART_UPLOAD_STATUS, MULTIPART_UPLOAD_STATUS as Status, statusIsError, statusIsUpload,
} from 'shared/utils/common/constants';
import { LoadingButton } from '@mui/lab';
import Clear from '@mui/icons-material/Clear';
import ReplayIcon from '@mui/icons-material/Replay';
import ProjectService from 'shared/services/Project.service';
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

  const color = status === 'DONE'
    ? 'lightGreen'
    : status === 'ABORTED' || (!submissionProgress && status === 'UPLOAD_PENDING')
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
                flexGrow: 1, display: 'flex', alignItems: 'center', justifyContent: 'space-between',
              }}
            >
              {statusIsUpload(submissionProgress.status)
                  && (
                  <>
                    <ProgressBar value={getSubmissionProgressPercentage(submissionProgress)} />
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
                {status === 'UPLOAD_PENDING'
                   && <Typography variant="caption">Upload failed or cancled</Typography>}
                {status === 'PROCESSING_PENDING'
                   && <ProcessingStatus />}
                {status === 'ABORTED'
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

          {/* <CustomButton
            type="primary"
            sx={{
              outerHeight: '70%',
            }}
          >
            Create Mapping
          </CustomButton> */}

        </Box>
      </Stack>
    </Card>
  );
}
