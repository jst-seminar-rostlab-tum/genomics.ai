import CenterFocusWeakIcon from '@mui/icons-material/CenterFocusWeak';
import ZoomInIcon from '@mui/icons-material/ZoomIn';
import ZoomOutIcon from '@mui/icons-material/ZoomOut';
import { Box, CircularProgress, IconButton } from '@mui/material';
import GeneMapperCategories from 'components/GeneMapper/Categories';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import Sidepanel from 'components/GeneMapper/Sidepanel';
import { resetZoom, zoomInN, zoomOutN } from 'components/Visualization/src/newZoom';
import { UmapVisualization2 } from 'components/Visualization/src/umapVisualization';
import { csv } from 'd3';
import React, {
  useEffect, useRef, useState,
} from 'react';
import { useParams } from 'react-router-dom';
import ProjectMock from 'shared/services/mock/projects';
import ProjectService from 'shared/services/Project.service';

/**
 * Shows the UMAP visualization for a given project.
 * @param projectId id of the project the result belongs to
 */
function GeneMapperResultView() {
  const [project, setProject] = useState(null);
  const umapContainer = useRef(null);
  const graphContainer = useRef(null);
  const [umap, setUmap] = useState(null);
  const [rendered, setRendered] = useState(false);
  const [umapSize, setUmapSize] = useState({
    width: 0,
    height: 0,
  });

  const { projectId } = useParams();

  useEffect(() => {
    ProjectService.getProject(projectId)
      .then((data) => setProject(data))
      .catch(() => {
        ProjectMock.getProject(1).then((data) => setProject(data));
      });
  }, [projectId]);

  useEffect(() => {
    if (project?.location) {
      csv(project.location).then((data) => {
        setUmap(new UmapVisualization2(umapContainer.current, data, graphContainer.current));
      });
    }
  }, [project]);

  useEffect(() => {
    if (umap && umapSize.width > 0 && umapSize.height > 0) {
      if (rendered) {
        umap.resize(umapSize.width, umapSize.height);
      } else {
        umap.render(umapSize.width, umapSize.height);
        setRendered(true);
      }
    }
  }, [umap, umapSize, rendered]);

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
                  hide={(category, value) => {
                    umap.after(category, value);
                  }}
                  show={() => {
                    umap.before();
                  }}
                />
              </Sidepanel>
              <Box
                sx={{
                  flexGrow: 1, display: 'flex', flexDirection: 'column', justifyContent: 'center', overflow: 'hidden',
                }}
              >
                <Box sx={{ display: 'flex', justifyContent: 'center' }}>
                  <IconButton onClick={() => zoomInN()}>
                    <ZoomInIcon />
                  </IconButton>
                  <IconButton onClick={() => zoomOutN()}>
                    <ZoomOutIcon />
                  </IconButton>
                  <IconButton onClick={() => resetZoom()}>
                    <CenterFocusWeakIcon />
                  </IconButton>
                </Box>
                <Box
                  sx={{
                    flexGrow: 1, display: 'flex', justifyContent: 'center', overflow: 'hidden',
                  }}
                  ref={umapContainer}
                />
              </Box>
              <Sidepanel title="Graphs" collapseToRight>
                <Box ref={graphContainer} />
              </Sidepanel>
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
