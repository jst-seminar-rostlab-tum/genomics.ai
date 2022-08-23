import { Box, CircularProgress } from '@mui/material';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import ResultVisualization from 'components/GeneMapper/ResultVisualization';
import React, {
  useEffect, useState,
} from 'react';
import { useParams } from 'react-router-dom';
import ProjectMock from 'shared/services/mock/projects';
import ProjectService from 'shared/services/Project.service';
import DemoService from 'shared/services/Demo.service';
import { PROJECT_STATUS } from 'shared/utils/common/constants';

/**
 * Shows the UMAP visualization for a given project.
 */
function GeneMapperResultView({ loggedIn = true }) {
  const [project, setProject] = useState(null);

  const { projectId } = useParams();

  useEffect(() => {
    // fetch the projects
    ProjectService.getProject(projectId)
      .then((data) => {
        let updatedData = data;
        console.log(JSON.stringify(data));
        DemoService.getDemos().then((demos) => {
          // loop over demo datasets and check whether the project is a demo
          demos.forEach((demo) => {
            const id = data.fileName.split('_')[2];
            if (id && id === demo._id) {
              // updating the demo dataset
              updatedData = {
                ...data,
                status: PROJECT_STATUS.DONE,
                location: demo.csvURL,
              };
            }
          });
          setProject(updatedData);
        });
      })
      .catch(() => {
        ProjectMock.getProject(projectId).then((data) => setProject(data));
      });
  }, [projectId]);

  return (
    <Box
      sx={{
        // 100vh - HeaderView header height - HeaderView content padding - Footer height
        height: 'calc(100vh - 75px - 24px - 40px)',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      {project
        ? (
          <>
            <GeneMapperResultHeader project={project} loggedIn={loggedIn} />
            <Box
              sx={{
                flexGrow: 1,
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'stretch',
                overflow: 'hidden',
              }}
            >
              <ResultVisualization dataUrl={project.location} loggedIn={loggedIn} />
            </Box>
          </>
        )
        : (
          <Box sx={{
            flexGrow: 1, display: 'flex', justifyContent: 'Center', alignItems: 'center',
          }}
          >
            <CircularProgress disableShrink />
          </Box>
        )}
    </Box>
  );
}

export default GeneMapperResultView;
