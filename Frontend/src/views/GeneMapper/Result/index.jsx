import { Box, CircularProgress } from '@mui/material';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import ResultVisualization from 'components/GeneMapper/ResultVisualization';
import React, {
  useEffect, useState,
} from 'react';
import { useParams } from 'react-router-dom';
import ProjectMock from 'shared/services/mock/projects';
import ProjectService from 'shared/services/Project.service';

/**
 * Shows the UMAP visualization for a given project.
 */
function GeneMapperResultView() {
  const [project, setProject] = useState(null);

  const { projectId } = useParams();

  useEffect(() => {
    ProjectService.getProject(projectId)
      .then((data) => setProject(data))
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
            <GeneMapperResultHeader project={project} />
            <Box
              sx={{
                flexGrow: 1,
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'stretch',
                overflow: 'hidden',
              }}
            >
              <ResultVisualization dataUrl={project.location} />
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
