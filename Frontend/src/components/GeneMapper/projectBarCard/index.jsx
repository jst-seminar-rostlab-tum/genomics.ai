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
import {
  ExpandLess, ExpandMore, InfoOutlined, StarBorder,
} from '@mui/icons-material';
import { Modal, ModalTitle } from 'components/Modal';
import TeamService from 'shared/services/Team.service';
import CustomButton from 'components/CustomButton';
import { TabCard } from '../TabCard';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import { GeneralCard } from 'components/GeneMapper/GeneralCardVariant';

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
  project, atlas, model, submissionProgress, cancelUpload, handleDelete, userTeams, addProjectToTeam,
}) {
  const history = useHistory();

  const color = project.status === PROJECT_STATUS.DONE
    ? 'lightGreen'
    : project.status === PROJECT_STATUS.ABORTED
    || (!submissionProgress && project.status === PROJECT_STATUS.UPLOAD_PENDING)
    || project.status === PROJECT_STATUS.PROCESSING_FAILED
    || submissionProgress?.status === MULTIPART_UPLOAD_STATUS.CANCELING
    || statusIsError(submissionProgress?.status)
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
    <>

      <GeneralCard>
        <CardActionArea disableTouchRipple sx={{ p: 0, height: '4em', borderRadius: '0.625rem' }}>
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
                   && <CanceldOrFailedStatus />}
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
          <ModalTitle>
            Why adding to a team?
          </ModalTitle>
          <Box sx={{ maxWidth: 320, overflow: 'auto' }}>
            <Typography>
              By Adding to a team, the visibility of
              a project will be inherited
              from the visibility of the team.
              If the assigned team is public,
              all teams will have access to the project.
              If the assigned team is private, the project will be accessible only to this team.
            </Typography>
          </Box>

        </Modal>
        <Box>
          {userTeams.map(
            (team) => (
              <TabCard
                data={{ name: team.name, visibility: team.visibility }}
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

    </>
  );
}
