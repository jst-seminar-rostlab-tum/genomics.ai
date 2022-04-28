import { Box, CircularProgress } from '@mui/material';
import GeneMapperCategories from 'components/GeneMapper/Categories';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import Sidepanel from 'components/GeneMapper/Sidepanel';
import { UmapVisualization2 } from 'components/Visualization/src/umapVisualization';
import { csv } from 'd3';
import React, {
  useEffect, useRef, useState,
} from 'react';
import getProject from 'shared/services/mock/projects';

/**
 * Shows the UMAP visualization for a given project.
 * @param projectId id of the project the result belongs to
 */
function GeneMapperResultView({ projectId }) {
  const [project, setProject] = useState(null);
  const umapContainer = useRef(null);
  const [umap, setUmap] = useState(null);
  const [umapSize, setUmapSize] = useState({
    width: 0,
    height: 0,
  });

  useEffect(() => {
    getProject(projectId)
      .then((data) => setProject(data));
  }, [projectId]);

  useEffect(() => {
    if (project?.resultURL) {
      csv(project.resultURL).then((data) => {
        setUmap(new UmapVisualization2(umapContainer.current, data));
      });
    }
  }, [project]);

  useEffect(() => {
    if (umap && umapSize.width > 0 && umapSize.height > 0) {
      umap.render(umapSize.width, umapSize.height);
    }
  }, [umap, umapSize]);

  useEffect(() => {
    if (umapContainer?.current) {
      const observer = new ResizeObserver((entries) => {
        const container = entries[0];
        setUmapSize({
          width: container.contentRect.width,
          height: container.contentRect.height,
        });
      });
      observer.observe(umapContainer.current);
      return () => observer.disconnect();
    }
    return () => {};
  }, [umapContainer.current]);

  return (
    <Box
      sx={{
        height: '100%',
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
                flexGrow: 1,
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'stretch',
                overflow: 'hidden',
              }}
            >
              <Sidepanel title="Categories">
                <GeneMapperCategories
                  categories={umap?.coloringModes}
                  setColorMode={(mode) => umap.setColorMode(mode)}
                />
              </Sidepanel>
              <Box
                sx={{
                  flexGrow: 1, display: 'flex', justifyContent: 'center', overflow: 'hidden',
                }}
                ref={umapContainer}
              />
              <Sidepanel title="Graphs" collapseToRight />
            </Box>
          </>
        )
        : (
          <Box sx={{
            flexGrow: 1, display: 'flex', justifyContent: 'Center', alignItems: 'center',
          }}
          >
            <CircularProgress />
          </Box>
        )}
    </Box>
  );
}

export default GeneMapperResultView;
