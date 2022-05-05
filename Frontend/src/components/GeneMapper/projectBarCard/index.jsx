import React, { useEffect, useState } from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import {
  Box, IconButton, LinearProgress, Stack, CardActionArea, FormControl, InputLabel, MenuItem, Select, Collapse, List, ListItemButton, ListItemIcon, ListItemText, Divider,
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
import { Modal, ModalTitle } from 'components/Modal';
import TeamService from 'shared/services/Team.service';
import CustomButton from 'components/CustomButton';
import { TabCard } from '../TabCard';

function ProcessingStatus() {
  return (
    <>
      <Box sx={{ pr: 2, flexGrow: 1 }}>
        <LinearProgress />
      </Box>
      <Typography variant="caption" noWrap sx={{ pr: 2 }}>Processing by scArches...</Typography>
    </>
  );
}

export default function ProjectBarCard({
  project, atlas, model, submissionProgress, setSubmissionProgress, handleDelete, userTeams, addProjectToTeam,
}) {
  const history = useHistory();

  const color = project.status === PROJECT_STATUS.DONE
    ? 'lightGreen'
    : project.status === PROJECT_STATUS.ABORTED
    || (!submissionProgress && project.status === PROJECT_STATUS.UPLOAD_PENDING)
    || project.status === PROJECT_STATUS.PROCESSING_FAILED
      ? 'red'
      : 'orange';

  const [projectTeam, setProjectTeam] = useState({});
  const [addTeam, setAddTeam] = useState(false);
  const [selectedTeam, setSelectedTeam] = useState('');
  const [open, setOpen] = useState(false);

  useEffect(() => {
    if (project.teamId) {
      TeamService.getTeam(project.teamId).then((team) => setProjectTeam(team));
    }
  }, [project.teamId]);

  const handleOpen = () => setAddTeam(true);
  const handleClose = () => setAddTeam(false);

  const handleClickCard = () => {
    setOpen(!open);
  };

  return (
    <>
      <Card sx={{
        marginTop: '5',
        marginBottom: '0.5em',
        borderStyle: 'solid',
        borderColor: '#C8C8C8',
        borderWidth: '0.1px',
      }}
      >

        <CardActionArea disableTouchRipple sx={{ p: 1 }}>
          <Stack direction="row" sx={{ justifyContent: 'space-between', alignItems: 'center' }}>

            <Stack
              direction="row"
            // spacing={4}
              sx={{ alignItems: 'center', flexGrow: 1 }}
              onClick={handleClickCard}
            >
              <Box sx={{ flexDirection: 'row', ml: 1, alignItems: 'center' }}>
                {open ? (
                  <ExpandLess sx={{
                    fontSize: 30,
                    transform: 'rotate(180deg)',
                  }}
                  />
                ) : (
                  <ExpandMore sx={{
                    fontSize: 30,
                    transform: 'rotate(-90deg)',
                  }}
                  />
                )}
              </Box>
              <CircleIcon sx={{
                fontSize: 30, color, mr: 1, ml: 1,
              }}
              />
              <Typography noWrap sx={{ width: '30%' }}>
                {project.name}
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
                    sx={{
                      flexGrow: 1, display: 'flex', alignItems: 'center',
                    }}
                  >
                    {project.status === PROJECT_STATUS.UPLOAD_PENDING
                   && (
                   <Typography variant="caption">
                     Upload failed or canceled
                   </Typography>
                   )}
                    {project.status === PROJECT_STATUS.PROCESSING_PENDING
                   && <ProcessingStatus />}
                    {(project.status === PROJECT_STATUS.ABORTED
                  || project.status === PROJECT_STATUS.PROCESSING_FAILED)
                   && <Typography variant="caption">Processing failed</Typography>}
                  </Box>
                )
                : null}
            </Stack>

            <Box sx={{
              p: 0.1, bgcolor: 'background.paper', borderRadius: 3, width: 'flex', mr: 1, display: 'flex', alignItems: 'center',
            }}
            >
              {projectTeam?.name
                ? (
                  <CustomButton type="tertiary" sx={{ mr: 1 }} onClick={() => history.push(`/sequencer/teams/${projectTeam._id || projectTeam.id}`)}>
                    <Typography>
                      {projectTeam.name}
                    </Typography>
                  </CustomButton>
                )
                : (
                  <Button
                    variant="outlined"
                    sx={{
                      borderRadius: 100,
                      mr: 1,
                    }}
                    style={{ textTransform: 'none' }}
                    onClick={handleOpen}
                  >
                    Add To Team
                  </Button>
                )}

              <Button
                variant="contained"
                color="secondary"
                sx={{
                  borderRadius: 100,
                  mr: 1,
                }}
                style={{ textTransform: 'none' }}
                onClick={() => history.push(`./genemapper/result/${project._id}`)}
                disabled={project.status !== 'DONE'}
              >
                See Results
              </Button>
              <IconButton onClick={() => handleDelete()}>
                <DeleteOutlineIcon color="error" />
              </IconButton>
            </Box>

          </Stack>
        </CardActionArea>

        <Collapse in={open} timeout="auto">
          <Divider variant="middle" />
          <Box sx={{ pl: 11.5, pb: 1, pt: 1 }}>
            <Typography>{`Atlas: ${atlas?.name}`}</Typography>
            <Typography>{`Model: ${model?.name}`}</Typography>
            <Typography>{`Dataset: ${project?.fileName}`}</Typography>
          </Box>
        </Collapse>

      </Card>
      <Modal
        isOpen={addTeam}
        setOpen={setAddTeam}
      >
        <ModalTitle>
          Select A Team
        </ModalTitle>
        <Box>
          {userTeams.map(
            (team) => (
              <TabCard
                data={{ name: team.name }}
                selected={team?._id === selectedTeam || team?.id === selectedTeam}
                handleOnClick={() => setSelectedTeam(team?._id || team?.id)}
              />
            ),
          )}
          <Box sx={{ display: 'flex', justifyContent: 'space-between', pt: 4 }}>
            <CustomButton
              type="tertiary"
              onClick={() => {
                setSelectedTeam('');
                setAddTeam(false);
              }}
            >
              Close

            </CustomButton>
            <CustomButton
              type="primary"
              onClick={() => {
                addProjectToTeam(selectedTeam);
                setAddTeam(false);
                setSelectedTeam('');
              }}
              disabled={selectedTeam === ''}
            >
              Confirm
            </CustomButton>
          </Box>
        </Box>
      </Modal>
    </>
  );
}
