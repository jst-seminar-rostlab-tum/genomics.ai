import * as React from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import {
  Box, IconButton, LinearProgress, Stack, CardActionArea, FormControl, InputLabel, MenuItem, Select, Collapse, List, ListItemButton, ListItemIcon, ListItemText,
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
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline';
import { ExpandLess, ExpandMore, StarBorder } from '@mui/icons-material';
import { Modal } from 'components/Modal';

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
  projectId, name, status, submissionProgress, setSubmissionProgress, handleDelete
}) {
  const history = useHistory();

  const color = status === PROJECT_STATUS.DONE
    ? 'lightGreen'
    : status === PROJECT_STATUS.ABORTED
    || (!submissionProgress && status === PROJECT_STATUS.UPLOAD_PENDING)
    || status === PROJECT_STATUS.PROCESSING_FAILED
      ? 'red'
      : 'orange';

  const [addTeam, setAddTeam] = React.useState(false);
  const [team, setTeam] = React.useState('');
  const handleOpen = () => setAddTeam(true);
  const handleClose = () => setAddTeam(false);
  const [open, setOpen] = React.useState(false);
  const handleChange = (event) => {
    setTeam(event.target.value);
    window.localStorage.setItem('Team', event.target.value);
  };

  const handleClickCard = () => {
    setOpen(!open);
  };

  React.useEffect(() => {
    if (window.localStorage.getItem('Team')) { setTeam(window.localStorage.getItem('Team')); }
  }, []);

  return (
    <Card sx={{
      marginTop: '5', marginBottom: '0.5em', borderStyle: 'solid', borderColor: '#C8C8C8', borderWidth: '0.1px',
    }}
    >

      <Stack direction="row" sx={{ justifyContent: 'space-between' }}>

        <CardActionArea disableTaouchRipple onClick={handleClickCard}>
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
            <Typography>
              {team}
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
                  sx={{ flexGrow: 1, display: 'flex' }}
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
            {open ? (
              <Box sx={{ flexDirection: 'row' }}>
                <ExpandLess sx={{
                  fontSize: 30, marginLeft: '3%',
                }}
                />
              </Box>
            ) : (
              <Box sx={{ flexDirection: 'row' }}>
                <ExpandMore sx={{
                  fontSize: 30, marginLeft: '3%',
                }}
                />
              </Box>
            )}

          </Stack>
          <Collapse in={open} timeout="auto" unmountOnExit collapsedSize="100px" sx={{}}>
            <Typography>More Infomation</Typography>
          </Collapse>
        </CardActionArea>

        <Box sx={{
          p: 0.1, bgcolor: 'background.paper', borderRadius: 3, width: 'flex', mr: 3, display: 'flex', alignItems: 'center',
        }}
        >
          <Button
            variant="outlined"
            size="small"
            sx={{
              borderRadius: 100,
              mr: 2,
              m: 1,
            }}
            style={{ textTransform: 'none' }}
            onClick={handleOpen}
          >
            Add To Team
          </Button>

          <Button
            variant="contained"
            color="secondary"
            sx={{
              borderRadius: 100,
              mr: 2,
            }}
            style={{ textTransform: 'none' }}
            onClick={() => history.push(`./genemapper/result/${projectId}`)}
            disabled={status !== 'DONE'}
          >
            See Results
          </Button>
          <IconButton>
            <DeleteOutlineIcon color="error" />
          </IconButton>
        </Box>

      </Stack>

      <Modal
        isOpen={addTeam}
        setOpen={setAddTeam}
        children={(
          <Box>
            <Typography variant="h6" sx={{ paddingBottom: 1 }}>
              Select A Team
            </Typography>
            <div>
              <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
                <InputLabel id="demo-simple-select-standard-label">Team</InputLabel>
                <Select
                  labelId="demo-simple-select-standard-label"
                  id="demo-simple-select-standard"
                  value={team}
                  onChange={handleChange}
                  label="Team"
                >
                  <MenuItem value="">
                    <em>None</em>
                  </MenuItem>
                  <MenuItem value="Team 1">Team 1</MenuItem>
                  <MenuItem value="Team 2">Team 2</MenuItem>
                  <MenuItem value="Team 3">Team 3</MenuItem>
                </Select>
              </FormControl>
            </div>
            <Button size="small" onClick={() => setAddTeam(false)}>Close</Button>
            <Button size="small" onClick={() => setAddTeam(false)}>Confirm</Button>
          </Box>
)}
      />

    </Card>
  );
}
