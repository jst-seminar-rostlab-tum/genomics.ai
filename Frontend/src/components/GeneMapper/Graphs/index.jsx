import React, { useEffect, useRef, useState } from 'react';
import { Box, IconButton, Modal } from '@mui/material';
import ZoomOutMapIcon from '@mui/icons-material/ZoomOutMap';
import GraphModal from '../GraphModal';

function GeneMapperGraphs({ graphs, drawGraph }) {
  const graphContainers = useRef({});
  const [open, setOpen] = useState('');

  useEffect(() => {
    Object.entries(graphContainers.current).forEach(([id, elem]) => {
      drawGraph(elem, id, 270, 270);
    });
  }, [graphContainers.current]);

  return (
    <>
      <Box sx={{ height: '100%' }}>
        {graphs.map((id) => (
          <Box
            key={id}
            ref={(elem) => { graphContainers.current[id] = elem; }}
            sx={{ position: 'relative' }}
          >
            <IconButton
              sx={{ position: 'absolute', top: 0, left: 0 }}
              onClick={() => {
                setOpen(id);
              }}
            >
              <ZoomOutMapIcon />
            </IconButton>
          </Box>
        ))}
      </Box>
      {open
        ? (
          <Modal open onClose={() => setOpen('')}>
            <Box bgcolor="red" height="500px" flex alignItems="center" justifyItems="center">
              <GraphModal draw={(container, w, h) => drawGraph(container, open, w, h)} />
            </Box>
          </Modal>
        )
        : null}
    </>
  );
}

export default GeneMapperGraphs;
