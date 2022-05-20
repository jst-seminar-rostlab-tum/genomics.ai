import React, { useEffect, useRef, useState } from 'react';
import { Box, IconButton } from '@mui/material';
import ZoomOutMapIcon from '@mui/icons-material/ZoomOutMap';
import { Modal } from 'components/Modal';

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

      <Modal isOpen={open} setOpen={() => setOpen('')}>
        <Box ref={(elem) => drawGraph(elem, open, 500, 500)} />
      </Modal>

    </>
  );
}

export default GeneMapperGraphs;
