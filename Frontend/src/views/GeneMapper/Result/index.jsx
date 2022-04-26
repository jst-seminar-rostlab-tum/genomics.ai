import { Box, CircularProgress } from '@mui/material';
import GeneMapperCategories from 'components/GeneMapper/Categories';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import Sidepanel from 'components/GeneMapper/Sidepanel';
import { csv } from 'd3';
import React, { useCallback, useEffect, useState } from 'react';
import getProject from 'shared/services/mock/projects';

const testCategories = {
  'Cell type': [
    {
      title: 'a',
      color: 'red',
    },
    {
      title: 'b',
      color: 'blue',
    },
    {
      title: 'c',
      color: 'green',
    },
  ],
  Batch: [
    {
      title: 'd',
      color: 'orange',
    },
    {
      title: 'e',
      color: 'yellow',
    },
    {
      title: 'f',
      color: 'lime',
    },
  ],
};

/**
 * Shows the UMAP visualization for a given project.
 * @param sidebarShown set true if sidebar is open, false otherwise
 */
function GeneMapperResultView({ sidebarShown, projectId }) {
  const paddingL = useCallback(() => (sidebarShown ? '350px' : '100px'), [sidebarShown]);

  const [project, setProject] = useState(null);

  useEffect(() => {
    getProject(projectId)
      .then((response) => response.json())
      .then((data) => setProject(data));
  }, [projectId]);

  useEffect(() => {
    if (project) { csv(project.resultURL).then((data) => console.log(data)); }
  }, [project]);

  return (
    <Box
      sx={{
        pl: paddingL,
        pr: '20px',
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      {project
        ? (
          <>
            <GeneMapperResultHeader projectName={project.name} />
            <Box
              sx={{
                display: 'flex', flexGrow: 1, justifyContent: 'space-between', alignItems: 'stretch',
              }}
            >
              <Sidepanel title="Categories">
                <GeneMapperCategories categories={testCategories} />
              </Sidepanel>
              <Box sx={{ flexGrow: 1 }} />
              <Sidepanel title="Graphs" collapseToRight />
            </Box>
          </>
        )
        : (
          <CircularProgress />
        )}
    </Box>
  );
}

export default GeneMapperResultView;
