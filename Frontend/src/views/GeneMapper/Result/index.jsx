import { Box, CircularProgress } from '@mui/material';
import GeneMapperCategories from 'components/GeneMapper/Categories';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import Sidepanel from 'components/GeneMapper/Sidepanel';
import { UmapVisualization2 } from 'components/Visualization/src/umapVisualization';
import { csv } from 'd3';
import React, {
  useCallback, useEffect, useRef, useState,
} from 'react';
import getProject from 'shared/services/mock/projects';

const testCategories = {
  cell_type: [
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
  batch: [
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
  const umapContainer = useRef(null);
  const [umap, setUmap] = useState(null);

  useEffect(() => {
    getProject(projectId)
      .then((response) => response.json())
      .then((data) => setProject(data));
  }, [projectId]);

  useEffect(() => {
    if (project?.resultURL) {
      console.log('new newmap');
      csv(project.resultURL).then((data) => {
        setUmap(new UmapVisualization2(umapContainer.current, data));
      });
    }
  }, [project]);

  useEffect(() => {
    if (umap) {
      console.log(umap, umapContainer.current);
      const height = umapContainer.current.clientHeight;
      const width = umapContainer.current.clientWidth;
      const dim = Math.min(height, width);
      umap.render(dim, dim);
    }
  }, [umap, umapContainer.current]);

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
                <GeneMapperCategories
                  categories={testCategories}
                  setColorMode={(mode) => umap.setColorMode(mode)}
                />
              </Sidepanel>
              <Box sx={{ flexGrow: 1 }} ref={umapContainer} />
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
