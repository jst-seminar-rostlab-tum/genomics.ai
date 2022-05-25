import React, { useEffect, useState } from 'react';
import { Modal, ModalTitle } from 'components/Modal';
import Button from 'components/CustomButton';
import {
  CircularProgress, Snackbar, Alert,
  DialogActions, DialogContent, DialogContentText,
} from '@mui/material';
import ProjectService from 'shared/services/Project.service';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import AtlasService from 'shared/services/Atlas.service';
import ModelService from 'shared/services/Model.service';
import TeamService from 'shared/services/Team.service';

function TeamProjectList({ team, isAdmin, updateTeam }) {
  const [projects, setProjects] = useState([]);
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [userTeams, setUserTeams] = useState([]);
  const [isLoading, setIsLoading] = useState(true);

  const [dialogOpen, setDialogOpen] = useState(false);
  const [snackbarOpen, setSnackbarOpen] = useState(false);
  const [removedProject, setRemovedProject] = useState([]);
  const [errorMessage, setErrorMessage] = useState('');

  const handleCloseDialog = () => setDialogOpen(false);

  function handleDelete(project) {
    if (isAdmin) {
      setDialogOpen(true);
      setRemovedProject(project);
    } else {
      setSnackbarOpen(true);
    }
  }

  async function remove() {
    setErrorMessage('');
    try {
      await TeamService.removeProjectFromTeam(team._id, removedProject._id);
      handleCloseDialog();
      updateTeam();
    } catch (e) {
      setErrorMessage(e.message);
    }
  }

  const handleSnackbarClose = (event, reason) => {
    if (reason === 'clickaway') {
      return;
    }
    setSnackbarOpen(false);
  };

  useEffect(() => {
    setIsLoading(true);
    AtlasService.getAtlases().then((data) => setAtlases(data));
    ModelService.getModels().then((data) => setModels(data));
    TeamService.getMyTeams().then((teams) => setUserTeams(teams));
    ProjectService.getTeamProjects(team._id).then((data) => setProjects(data));
    setIsLoading(false);
  }, []);

  if (isLoading) {
    return <CircularProgress />;
  }

  if (projects.length === 0) {
    return (<span>No projects.</span>);
  }

  return (
    <div>
      {projects
        .map((project) => (
          <ProjectBarCard
            key={project._id}
            project={project}
            atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
            model={models.find((model) => String(model._id) === String(project.modelId))}
            userTeams={userTeams}
            handleDelete={() => handleDelete(project)}
          />
        ))}
      <Modal
        isOpen={dialogOpen}
        setOpen={(o) => !o && handleCloseDialog()}
      >
        <ModalTitle>
          Remove Project
        </ModalTitle>
        <DialogContent>
          <DialogContentText>
            {`Do you really want to remove the project ${removedProject.name} from the team ${team.name}?`}
          </DialogContentText>
          {
            errorMessage && (
              <DialogContentText id="alert-dialog-description" color="error">
                {errorMessage}
              </DialogContentText>
            )
          }
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDialog} type="tertiary">Cancel</Button>
          <Button onClick={() => remove()} type="critical" autoFocus>
            Remove
          </Button>
        </DialogActions>
      </Modal>
      <Snackbar
        open={snackbarOpen}
        autoHideDuration={3000}
        onClose={handleSnackbarClose}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'center',
        }}
      >
        <Alert onClose={handleSnackbarClose} severity="error" sx={{ width: '100%' }}>
          Only an admin may remove projects from a team!
        </Alert>
      </Snackbar>
    </div>
  );
}

export default TeamProjectList;
