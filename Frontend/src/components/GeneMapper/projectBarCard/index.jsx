import React, { useEffect, useState } from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import {
  Box, IconButton, LinearProgress, Stack, CardActionArea, FormControl, InputLabel, MenuItem, Select, Collapse, List, ListItemButton, ListItemIcon, ListItemText, Divider, Grid, Alert,
} from '@mui/material';
import CircleIcon from '@mui/icons-material/Circle';
import DownloadIcon from '@mui/icons-material/Download';
import { useHistory } from 'react-router-dom';
import { getSubmissionProgressPercentage } from 'shared/services/UploadLogic';
import {
  MULTIPART_UPLOAD_STATUS, MULTIPART_UPLOAD_STATUS as Status, statusIsError, statusIsUpload, PROJECT_STATUS,
} from 'shared/utils/common/constants';
import Clear from '@mui/icons-material/Clear';
import ReplayIcon from '@mui/icons-material/Replay';
import ProgressBar from 'components/ProgressBar';
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline';
import {
  ExpandLess, ExpandMore, InfoOutlined, StarBorder,
} from '@mui/icons-material';
import { Modal, ModalTitle } from 'components/Modal';
import TeamService from 'shared/services/Team.service';
import CustomButton from 'components/CustomButton';
import { TabCard } from '../TabCard';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import { GeneralCard } from 'components/Cards/GeneralCard';
import ProjectInfo from '../ProjectInfo';
import { initSubmissionProgress, useSubmissionProgress } from 'shared/context/submissionProgressContext';

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

function CanceldOrFailedStatus() {
  return (
    <Typography variant="caption">
      Upload failed or canceled
    </Typography>
  );
}

export default function ProjectBarCard({
  project, atlas, model, userTeams, handleDelete, deleted,
}) {
  const history = useHistory();
  const [submissionProgresses, setSubmissionProgresses] = useSubmissionProgress();

  const submissionProgress = submissionProgresses[project._id];

  const cancelUpload = () => {
    setSubmissionProgresses((prev) => ({
      ...prev,
      [project._id]: {
        ...(prev[project._id] ?? initSubmissionProgress(project.uploadId)),
        status: MULTIPART_UPLOAD_STATUS.CANCELING,
      },
    }));
    localStorage.setItem(`cancelUpload_${project.uploadId}`, '1');
  };

  const addProjectToTeam = async (teamId) => {
    TeamService.addProject(teamId, project._id);
  };

  const color = project.status === PROJECT_STATUS.DONE
    ? 'lightGreen'
    : project.status === PROJECT_STATUS.ABORTED
    || (!submissionProgress && project.status === PROJECT_STATUS.UPLOAD_PENDING)
    || project.status === PROJECT_STATUS.PROCESSING_FAILED
    || submissionProgress?.status === MULTIPART_UPLOAD_STATUS.CANCELING
    || (submissionProgress && statusIsError(submissionProgress.status))
      ? 'red'
      : 'orange';

  const [projectTeam, setProjectTeam] = useState({});
  const [addTeam, setAddTeam] = useState(false);
  const [selectedTeam, setSelectedTeam] = useState('');
  const [open, setOpen] = useState(false);
  const [showInfo, setShowInfo] = useState(false);

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
    <Box sx={{ mb: 2 }}>
      <GeneralCard padding={0}>
        <CardActionArea
          disableTouchRipple
          sx={{
            p: 2, borderRadius: 'inherit',
          }}
        >
          <Grid container direction="row" justifyContent="space-between" sx={{ width: '100%', justifyContent: 'space-between', alignItems: 'center' }}>
            <Grid
              item
              xs={12}
              md={8}
              container
              direction="row"
              alignItems="center"
              onClick={handleClickCard}
            >
              <Grid container item xs={4} alignItems="center">
                <Box sx={{ flexDirection: 'row', alignItems: 'center' }}>
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
                  fontSize: 30, color, mr: 1,
                }}
                />
                <Typography noWrap>
                  {project.name}
                </Typography>
              </Grid>
              <Grid item xs={8}>
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
                      onClick={(e) => {
                        e.stopPropagation();
                        cancelUpload();
                      }}
                    >
                      <Clear color="error" />
                    </IconButton>
                  </>
                  )}
                    {statusIsError(submissionProgress.status)
                  && <CanceldOrFailedStatus />}
                    {submissionProgress.status === MULTIPART_UPLOAD_STATUS.CANCELING
                  && <CanceldOrFailedStatus />}
                    {submissionProgress.status === MULTIPART_UPLOAD_STATUS.COMPLETE
                   && project.status !== PROJECT_STATUS.DONE
                   && project.status !== PROJECT_STATUS.ABORTED
                   && project.status !== PROJECT_STATUS.PROCESSING_FAILED
                   && <ProcessingStatus />}
                    {(project.status === PROJECT_STATUS.ABORTED
                  || project.status === PROJECT_STATUS.PROCESSING_FAILED)
                   && <Typography variant="caption">Processing failed</Typography>}
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
                   && <CanceldOrFailedStatus />}
                      {project.status === PROJECT_STATUS.PROCESSING_PENDING
                   && <ProcessingStatus />}
                      {(project.status === PROJECT_STATUS.ABORTED
                  || project.status === PROJECT_STATUS.PROCESSING_FAILED)
                   && <Typography variant="caption">Processing failed</Typography>}
                    </Box>
                  )
                  : null}
              </Grid>
            </Grid>
            <Grid container item md={4} justifyContent="flex-end">
              <Box sx={{
                p: 0.1, bgcolor: 'background.paper', borderRadius: 3, width: 'flex', display: 'flex', alignItems: 'center',
              }}
              >
                {!deleted
                && (
                <>
                  {projectTeam?.title
                    ? (
                      <CustomButton type="tertiary" sx={{ mr: 1 }} onClick={() => history.push(`/sequencer/teams/${projectTeam._id || projectTeam.id}`)}>
                        <Typography>
                          {projectTeam.title}
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
                    onClick={() => history.push(`/sequencer/genemapper/result/${project._id}`)}
                    disabled={project.status !== 'DONE'}
                  >
                    See Results
                  </Button>
                  <IconButton
                    href={project.location}
                    download={`${project.name}.tsv`}
                    disabled={project.status !== 'DONE'}
                  >
                    <DownloadIcon />
                  </IconButton>
                </>
                )}
                <IconButton onClick={() => handleDelete()}>
                  {deleted
                    ? <ReplayIcon />
                    : <DeleteOutlineIcon color="error" />}
                </IconButton>
              </Box>
            </Grid>
          </Grid>
        </CardActionArea>

        <Collapse in={open} timeout="auto">
          <Divider variant="middle" />
          <Box sx={{ pl: 11.5, pb: 1, pt: 1 }}>
            <ProjectInfo project={project} atlas={atlas} model={model} />
          </Box>
        </Collapse>

      </GeneralCard>
      <Modal
        isOpen={addTeam}
        setOpen={setAddTeam}
        sx={{ position: 'fixed', top: '20%' }}
      >
        <ModalTitle>
          <Stack direction="row" alignItems="center">
            Select a Team
            <IconButton size="small" onClick={() => setShowInfo(true)}>
              <InfoOutlined fontSize="small" />
            </IconButton>
          </Stack>
        </ModalTitle>
        <Modal
          isOpen={showInfo}
          setOpen={setShowInfo}
          sx={{ position: 'fixed', top: '20%' }}
        >
          <Box sx={{
            display: 'flex',
            flexDirection: 'column',
            alignItems: 'flex-start',
            justifyContent: 'space-between',
            ml: 3,
            mr: 3,
          }}
          >
            <Box sx={{
              display: 'flex', flexDirection: 'row', width: '100%', justifyContent: 'center',
            }}
            >
              <ModalTitle>
                Why Adding to a Team?
              </ModalTitle>
            </Box>

            <Box>
              <Typography sx={{ width: '100%', maxWidth: '800px' }}>
                By Adding to a team, the visibility of
                a project will be inherited
                from the visibility of the team.
                If the assigned team is public,
                all teams will have access to the project.
                If the assigned team is private, the project will be
                accessible only to this team.
              </Typography>
            </Box>
          </Box>

        </Modal>
        { userTeams?.length === 0
        && (
        <Alert severity="info">
          You have no existing teams. Please add your teams in community.
        </Alert>
        )}
        <Box>
          {userTeams.map(
             (team) => (
               <TabCard
                 key={team._id}
                 data={{ name: team.title, visibility: team.visibility.toLowerCase() }}
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
                TeamService.getTeam(selectedTeam).then((team) => setProjectTeam(team));
              }}
              disabled={selectedTeam === ''}
            >
              Confirm
            </CustomButton>
          </Box>
        </Box>
      </Modal>

    </Box>
  );
}
