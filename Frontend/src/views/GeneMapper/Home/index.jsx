import React, { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import {
  Typography, createTheme, ThemeProvider, Stack, TextField, Alert, CircularProgress,
} from '@mui/material';
import PlusIcon from 'components/general/PlusIcon';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import SearchIcon from '@mui/icons-material/Search';
import ProjectService from 'shared/services/Project.service';
import { initSubmissionProgress, useSubmissionProgress } from 'shared/context/submissionProgressContext';
import { MULTIPART_UPLOAD_STATUS, PROJECTS_UPDATE_INTERVAL, statusIsError } from 'shared/utils/common/constants';
import ProjectMock from 'shared/services/mock/projects';
import AtlasService from 'shared/services/Atlas.service';
import ModelService from 'shared/services/Model.service';
import TeamService from 'shared/services/Team.service';
import { useHistory } from 'react-router-dom';

const theme = createTheme({
  palette: {
    primary: {
      main: '#000000',
    },
    secondary: {
      main: '#5676E4',
    },
  },
});

function GeneMapperHome() {
  const [projects, setProjects] = useState(null);
  const [deletedProjects, setDeletedProjects] = useState([]);
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [userTeams, setUserTeams] = useState([]);

  const [findString, setFindString] = useState('');

  const history = useHistory();

  const handleDeleteProject = (id) => {
    // setProjects(projects.filter((object) => object._id != id));
    // const deleted = window.localStorage.getItem('DeletedProjects') ?? '';
    // console.log(deleted);

    // window.localStorage.setItem('DeletedProjects', `${deleted},${id}`);
    // ProjectMock.deleteProject(id);
    ProjectService.deleteProject(id).then(() => {
      ProjectService.getOwnProjects().then((data) => setProjects(data));
      ProjectService.getDeletedProjects().then((data) => setDeletedProjects(data));
    });
  };

  const handleRestoreProject = (id) => {
    ProjectService.restoreProject(id).then(() => {
      ProjectService.getOwnProjects().then((data) => setProjects(data));
      ProjectService.getDeletedProjects().then((data) => setDeletedProjects(data));
    });
  };

  useEffect(() => {
    ProjectService.getOwnProjects().then((data) => setProjects(data));
    const timer = setInterval(() => {
      ProjectService.getOwnProjects().then((data) => setProjects(data));
    }, PROJECTS_UPDATE_INTERVAL);

    return () => {
      clearInterval(timer);
    };
  }, []);

  useEffect(() => {
    AtlasService.getAtlases().then((data) => setAtlases(data));
    ModelService.getModels().then((data) => setModels(data));
    TeamService.getMyTeams().then((teams) => setUserTeams(teams));
    ProjectService.getDeletedProjects().then((data) => setDeletedProjects(data));
  }, []);

  return (
    <div>
      <ThemeProvider theme={theme}>
        <Box
          sx={{
            display: 'flex',
            flexDirection: 'row',
            justifyContent: 'space-between',
            pt: 2,
            pb: 3,
          }}
        >
          <Stack direction="row" className="stack" alignItems="Center">

            <Typography variant="h5" sx={{ pr: 1 }}>Your Mappings</Typography>
            <PlusIcon onClick={() => { history.push('genemapper/create'); }} />
          </Stack>
          <TextField
            id="outlined-basic"
            sx={{ width: '32.7ch' }}
            label={(
              <Stack direction="row">
                <SearchIcon />
                Find a Mapping
              </Stack>
            )}
            variant="outlined"
            size="small"
            value={findString}
            onChange={(e) => setFindString(e.target.value)}
          />
        </Box>
        {projects === null
        && <CircularProgress />}
        { projects === []
        && (
        <Alert severity="info">
          You have not created any mappings yet. Create one by clicking the Plus-Icon or learn more about ScArches by clicking the Help-Icon next to the title.
        </Alert>
        )}
        {projects
        && (
        <div>
          {projects.map((project) => (
            <ProjectBarCard
              key={project._id}
              project={project}
              atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
              model={models.find((model) => String(model._id) === String(project.modelId))}
              userTeams={userTeams}
              handleDelete={() => handleDeleteProject(project._id)}
            />

          ))}
        </div>
        )}
        {deletedProjects.length > 0
        && (
        <Box>
          <Typography variant="h6" sx={{ mt: 4, mb: 2 }}>Deleted Projects</Typography>
          {deletedProjects.map((project) => (
            <ProjectBarCard
              key={project._id}
              project={project}
              atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
              model={models.find((model) => String(model._id) === String(project.modelId))}
              userTeams={userTeams}
              handleDelete={() => handleRestoreProject(project._id)}
            />
          ))}
        </Box>
        )}
      </ThemeProvider>
    </div>
  );
}

export default GeneMapperHome;
