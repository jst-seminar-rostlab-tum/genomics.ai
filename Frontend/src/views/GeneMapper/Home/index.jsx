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
import {
  MULTIPART_UPLOAD_STATUS, PROJECTS_UPDATE_INTERVAL, PROJECT_STATUS, statusIsError,
} from 'shared/utils/common/constants';
import ProjectMock from 'shared/services/mock/projects';
import AtlasService from 'shared/services/Atlas.service';
import ModelService from 'shared/services/Model.service';
import DemoService from 'shared/services/Demo.service';
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

// The loggedIn parameter is set to true when the user is logged in
function GeneMapperHome({ loggedIn }) {
  const [projects, setProjects] = useState(null);
  const [deletedProjects, setDeletedProjects] = useState([]);
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);
  const [userTeams, setUserTeams] = useState([]);
  const [demoDatasets, setDemoDatasets] = useState([]);

  const [findString, setFindString] = useState('');

  const history = useHistory();

  /**
   * Function to find matching demos in the data variable.
   * If demos are found, the information stored about them in the "data" variable is updated.
   * @param data: The data variable containing the projects to search for demos in.
   * @param demos: The array of available demos
   */
  const findDemos = (data, demos) => {
    // loop over projects and search for demo datasets and set their attributes
    for (let i = 0; i < data.length; i++) {
      for (let j = 0; j < demos.length; j++) {
        // If the project is a demo, the id is stored in the file name.
        // file_name: <atlas>_<model>_<demo_id>
        const id = data[i].fileName.split('_')[2];
        // update the stored information about the demo datasets
        if (id && id === demos[j]._id) {
          // updating the demo dataset
          data[i] = {
            ...data[i],
            status: PROJECT_STATUS.DONE,
            location: demos[j].csvURL,
          };
        }
      }
    }
  };

  // function to handle deletion of project
  const handleDeleteProject = (id) => {
    DemoService.getDemos().then((demos) => {
      if (loggedIn) {
        ProjectService.deleteProject(id).then(() => {
          ProjectService.getOwnProjects().then((data) => {
            // search for demos and set the information stored about them
            findDemos(data, demos);
            setProjects(data);
          });
          ProjectService.getDeletedProjects().then((data) => {
            // search for demos and set the information stored about them
            findDemos(data, demos);
            setDeletedProjects(data);
          });
        });
      } else {
        // todo: delete the project from the cache and from the local storage
      }
    });
  };

  // function to handle restoration of project
  const handleRestoreProject = (id) => {
    DemoService.getDemos().then((demos) => {
      if (loggedIn) {
        ProjectService.restoreProject(id).then(() => {
          ProjectService.getOwnProjects().then((data) => {
            // search for demos and set the information stored about them
            findDemos(data, demos);
            setProjects(data);
          });
          ProjectService.getDeletedProjects().then((data) => {
            // search for demos and set the information stored about them
            findDemos(data, demos);
            setDeletedProjects(data);
          });
        });
      } else {
        //todo: fetch from the cache
      }
    });
  };

  // in the cache, store an array of projects using the following object structure
  // {_id, owner, name, modelId, atlasId, fileName, fileSize, status, resultSize, _v, uploadId, location}
  // the location is going to be the cache.
  // TODO: add the above to the getProjects function
  // link to get the cache:
  // https://www.geeksforgeeks.org/how-to-get-complete-cache-data-in-reactjs/

  // Function to get projects and update the necessary info about the demo datasets
  const getProjects = () => {
    // fetch demos
    DemoService.getDemos().then((demos) => {
      if (loggedIn) {
        ProjectService.getOwnProjects().then((data) => {
          // search for demos and set the information stored about them
          findDemos(data, demos);
          setDemoDatasets(demos);
          setProjects(data);
        });
      } else { // if not logged in, check the local storage and cache for the projects
        // todo: add the cache check for the projects
      }
    });
  };

  useEffect(() => {
    getProjects();
    const timer = setInterval(() => getProjects(), PROJECTS_UPDATE_INTERVAL);

    return () => {
      clearInterval(timer);
    };
  }, []);

  useEffect(() => {
    AtlasService.getAtlases().then((data) => setAtlases(data));
    ModelService.getModels().then((data) => setModels(data));
    if (loggedIn) {
      TeamService.getMyTeams().then((teams) => setUserTeams(teams));
      ProjectService.getDeletedProjects().then((data) => setDeletedProjects(data));
    }
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
          && (
            <Box sx={{ textAlign: 'center' }}>
              <CircularProgress />
            </Box>
          )}
        {/* todo: in the future, check the cache for the projects
        instead of disabling if not logged in, right now,
        it is disabled in order to not get null exceptions */}
        {loggedIn && projects?.length === 0
          && (
            <Alert severity="info">
              You have not created any mappings yet.
              Create one by clicking the Plus-Icon
              or learn more about ScArches by clicking the Help-Icon next to the title.
            </Alert>
          )}
        {/* todo: in the future, check the cache for the projects
        instead of disabling if not logged in, right now,
        it is disabled in order to not get null exceptions */}
        {loggedIn && projects
          && (
            <div>
              {projects.filter((project) => (
                (findString === '' || project.name.toLowerCase().includes(findString.toLowerCase())))).map((project) => (
                  <ProjectBarCard
                    key={project._id}
                    project={project}
                    atlas={atlases.find((atlas) => String(atlas._id) === String(project.atlasId))}
                    model={models.find((model) => String(model._id) === String(project.modelId))}
                    userTeams={userTeams}
                    handleDelete={() => handleDeleteProject(project._id)}
                    loggedIn={loggedIn}
                  />
              ))}
            </div>
          )}
        {loggedIn && deletedProjects.length > 0
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
                  loggedIn={loggedIn}
                  handleDelete={() => handleRestoreProject(project._id)}
                  deleted
                />
              ))}
            </Box>
          )}
      </ThemeProvider>
    </div>
  );
}

export default GeneMapperHome;
