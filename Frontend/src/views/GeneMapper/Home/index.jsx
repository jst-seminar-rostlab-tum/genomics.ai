import React, { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import {
  Typography, createTheme, ThemeProvider, Stack, TextField,
} from '@mui/material';
import PlusIcon from 'components/general/PlusIcon';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import SearchIcon from '@mui/icons-material/Search';
import ProjectService from 'shared/services/Project.service';
import { useSubmissionProgress } from 'shared/context/submissionProgressContext';
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
  const [projects, setProjects] = useState([]);
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [userTeams, setUserTeams] = useState([]);

  const [findString, setFindString] = useState('');
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();

  const history = useHistory();

  const handleDeleteItem = (id) => {
    setProjects(projects.filter((object) => object._id != id));
    const deleted = window.localStorage.getItem('DeletedProjects') ?? '';
    console.log(deleted);

    window.localStorage.setItem('DeletedProjects', `${deleted},${id}`);
    // ProjectMock.deleteProject(id);
  };

  const addProjectToTeam = async (teamId, projectId) => {
    const projectTeams = JSON.parse(window.localStorage.getItem('projectTeams')) || {};
    window.localStorage.setItem('projectTeams', JSON.stringify({ ...projectTeams, [projectId]: teamId }));
  };

  const teamOfProject = (projectId) => {
    const projectTeams = JSON.parse(window.localStorage.getItem('projectTeams')) || {};
    return projectTeams[projectId];
  };

  useEffect(() => {
    ProjectService.getOwnProjects().then((data) => setProjects(data));
    const timer = setInterval(() => {
      ProjectService.getOwnProjects().then((data) => setProjects(data));
      if (submissionProgress.status === MULTIPART_UPLOAD_STATUS.COMPLETE
        || submissionProgress.status === MULTIPART_UPLOAD_STATUS.CANCELING
        || statusIsError(submissionProgress.status)) {
        setSubmissionProgress({
          status: MULTIPART_UPLOAD_STATUS.IDLE,
          uploadId: '',
          chunks: 0,
          uploaded: 0,
          remaining: [],
          uploadedParts: [],
        });
      }
    }, PROJECTS_UPDATE_INTERVAL);

    return () => {
      clearInterval(timer);
    };
  }, [submissionProgress.status]);

  useEffect(() => {
    AtlasService.getAtlases().then((data) => setAtlases(data));
    ModelService.getModels().then((data) => setModels(data));
    TeamService.getMyTeams().then((teams) => setUserTeams(teams));
  }, []);

  return (
    <div>
      <ThemeProvider theme={theme}>
        {/* {Object.entries(submissionProgress)
          .map(([key, value]) => <Typography>{`${key}: ${value}` }</Typography>)} */}
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
        <div>
          {projects
            .filter((project) => (
              (findString === '' || project.name.toLowerCase().includes(findString.toLowerCase())))
              && !(window.localStorage.getItem('DeletedProjects') ?? []).includes(project._id))
            .map((project) => {
              const projectTeamId = teamOfProject(project._id);
              return (
                <ProjectBarCard
                  key={project._id}
                  project={projectTeamId ? { ...project, teamId: projectTeamId } : project}
                  atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
                  model={models.find((model) => String(model._id) === String(project.modelId))}
                  userTeams={userTeams}
                  addProjectToTeam={(teamId) => addProjectToTeam(teamId, project._id)}
                  handleDelete={() => handleDeleteItem(project._id)}
                  submissionProgress={submissionProgress.uploadId === project.uploadId
                    ? submissionProgress : null}
                  setSubmissionProgress={submissionProgress.uploadId === project.uploadId
                    ? setSubmissionProgress : () => { }}
                />

              );
            })}
        </div>

      </ThemeProvider>
    </div>
  );
}

export default GeneMapperHome;
