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
import CancelIcon from '@mui/icons-material/Cancel';
import { MULTIPART_UPLOAD_STATUS as Status, statusIsError, statusIsUpload } from 'shared/utils/common/constants';
import { LoadingButton } from '@mui/lab';
import Clear from '@mui/icons-material/Clear';
import ReplayIcon from '@mui/icons-material/Replay';
import ProjectService from 'shared/services/Project.service';
import ProgressBar from 'components/ProgressBar';

export default function ProjectBarCard({
  projectId, name, status, submissionProgress, setSubmissionProgress,
}) {
  const history = useHistory();
  const [color, setColor] = React.useState(status === 'DONE' ? 'lightGreen' : status === 'IN PROGRESS' ? 'orange' : 'red');
  const [typographyColor, setTypographyColor] = React.useState(status === 'UPLOAD FAILED' ? 'red' : 'black');

  return (
    <Card sx={{
      marginTop: '5', marginBottom: '0.5em', borderStyle: 'solid', borderColor: '#C8C8C8', borderWidth: '0.1px',
    }}
    >
      <Stack direction="row" sx={{ justifyContent: 'space-between', height: '56px' }}>
        <Stack
          direction="row"
          spacing={4}
          sx={{ width: '50%', alignItems: 'center' }}
        >

          <CircleIcon sx={{
            fontSize: 30, marginLeft: '3%', color,
          }}
          />
          <Typography noWrap sx={{ width: '30%' }}>
            {name}
          </Typography>
          {/* <Typography sx={{ color: typographyColor }}>
            {status}
          </Typography> */}
          {submissionProgress ? (
            <Box
              sx={{ flexGrow: 1, display: 'flex', alignItems: 'center' }}
            >
              {statusIsUpload(submissionProgress.status) ? (
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
              )
                : statusIsError(submissionProgress.status)
                  ? (
                    <>
                      <Typography>{submissionProgress.status}</Typography>
                      <IconButton onClick={() => {}}>
                        <ReplayIcon />
                      </IconButton>
                    </>
                  ) : null}
            </Box>
          ) : null}
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
