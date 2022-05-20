import CenterFocusWeakIcon from '@mui/icons-material/CenterFocusWeak';
import ZoomInIcon from '@mui/icons-material/ZoomIn';
import ZoomOutIcon from '@mui/icons-material/ZoomOut';
import { Box, CircularProgress, IconButton } from '@mui/material';
import { resetZoom, zoomInN, zoomOutN } from 'components/Visualization/src/newZoom';
import { UmapVisualization2 } from 'components/Visualization/src/umapVisualization';
import { csv } from 'd3';
import React, { useEffect, useRef, useState } from 'react';
import GeneMapperCategories from '../Categories';
import Sidepanel from '../Sidepanel';

/**
 *
 * @param dataUrl download link for the csv data
 * @param onlyUmap set to true if only the UMAP should be shown
 * @dataUrl url to download the data from
 */
function ResultVisualization({ dataUrl, onlyUmap }) {
  const umapContainer = useRef(null);
  const graphContainer = useRef(null);

  const [umap, setUmap] = useState(null);
  const [umapSize, setUmapSize] = useState({
    width: 0,
    height: 0,
  });

  useEffect(() => {
    csv(dataUrl).then((data) => {
      setUmap(new UmapVisualization2(umapContainer.current, data, graphContainer.current));
    });
  }, [dataUrl]);

  useEffect(() => {
    if (umap && umapSize.width > 0 && umapSize.height > 0) {
      umap.resize(umapSize.width, umapSize.height);
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
    <Box sx={{ display: 'flex', width: '100%', height: '100%' }}>
      {!umap ? (
        <Box sx={{
          position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%, -50%)',
        }}
        >
          <CircularProgress disableShrink />
        </Box>
      ) : null}
      <Box sx={{ display: onlyUmap ? 'none' : 'block' }}>
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
      </Box>
      <Box
        sx={{
          flexGrow: 1,
          display: 'flex',
          flexDirection: 'column',
          justifyContent: 'center',
          overflow: 'hidden',
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
      <Box sx={{ display: onlyUmap ? 'none' : 'block' }}>
        <Sidepanel title="Graphs" collapseToRight>
          <Box ref={graphContainer} />
        </Sidepanel>
      </Box>
    </Box>
  );
}

export default ResultVisualization;
