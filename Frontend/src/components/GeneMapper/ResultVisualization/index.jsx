import CenterFocusWeakIcon from '@mui/icons-material/CenterFocusWeak';
import ZoomInIcon from '@mui/icons-material/ZoomIn';
import ZoomOutIcon from '@mui/icons-material/ZoomOut';
import PublicIcon from '@mui/icons-material/Public';
import BiotechIcon from '@mui/icons-material/Biotech';
import BlurOnIcon from '@mui/icons-material/BlurOn';
import {
  Box, CircularProgress, Divider, IconButton, Tooltip,
} from '@mui/material';
import { resetZoom, zoomInN, zoomOutN } from 'components/Visualization/src/newZoom';
import { UmapVisualization2 } from 'components/Visualization/src/umapVisualization';
import { csv } from 'd3';
import React, { useEffect, useRef, useState } from 'react';
import GeneMapperCategories from '../Categories';
import Sidepanel from '../Sidepanel';
import { colors } from 'shared/theme/colors';
import GeneMapperGraphs from '../Graphs';

const activatedColor = colors.primary['400'];
const deactivatedColor = colors.primary['200'];

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

  const [showQuery, setShowQuery] = useState(true);
  const [showReference, setShowReference] = useState(true);
  const [showPredictions, setShowPredictions] = useState(true);

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
            show={(category, value) => {
              umap.before(category, value);
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
          {umap
          && Object.keys(umap.coloringModes).includes('predicted')
          && (
            <>
              <Box>
                <Tooltip title={`${showPredictions ? 'Hide' : 'Show'} predicted data`}>
                  <IconButton onClick={() => {
                    setShowPredictions(!showPredictions);
                    if (showPredictions) {
                      umap.predictedCellsTransparent();
                    } else {
                      umap.predictedCellsVisible();
                    }
                  }}
                  >
                    <BlurOnIcon
                      sx={{ color: showPredictions ? activatedColor : deactivatedColor }}
                    />
                  </IconButton>
                </Tooltip>
              </Box>
              <Divider orientation="vertical" flexItem variant="middle" />
            </>
          )}
          {
            umap
            && Object.keys(umap.coloringModes).includes('type')
            && (
              <>
                <Box>
                  <Tooltip title={`${showReference ? 'Hide' : 'Show'} reference data`}>
                    <IconButton onClick={() => {
                      setShowReference(!showReference);
                      if (showReference) {
                        setShowQuery(true);
                        umap.showQuery();
                        umap.hideReference();
                      } else {
                        umap.showReference();
                      }
                    }}
                    >
                      <PublicIcon
                        sx={{ color: showReference ? activatedColor : deactivatedColor }}
                      />
                    </IconButton>
                  </Tooltip>
                  <Tooltip title={`${showQuery ? 'Hide' : 'Show'} query data`}>
                    <IconButton onClick={() => {
                      setShowQuery(!showQuery);
                      if (showQuery) {
                        setShowReference(true);
                        umap.showReference();
                        umap.hideQuery();
                      } else {
                        umap.showQuery();
                      }
                    }}
                    >
                      <BiotechIcon sx={{ color: showQuery ? activatedColor : deactivatedColor }} />
                    </IconButton>
                  </Tooltip>
                </Box>
                <Divider orientation="vertical" flexItem variant="middle" />
              </>
            )
          }
          <Box>
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
          {umap
          && (
          <GeneMapperGraphs
            graphs={umap.getAvailableGraphs()}
            drawGraph={(container, graph, width, height) => {
              umap.drawGraph(container, graph, width, height);
            }}
          />
          )}
        </Sidepanel>
      </Box>
    </Box>
  );
}

export default ResultVisualization;
