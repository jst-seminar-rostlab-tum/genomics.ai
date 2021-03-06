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
              size="small"
              sx={{ position: 'absolute', left: 50 }}
              onClick={() => {
                setOpen(id);
              }}
            >
              <ZoomOutMapIcon />
            </IconButton>
          </Box>
        ))}
      </Box>

      <Modal isOpen={open} setOpen={() => setOpen('')} maxWidth="false">
        <Box ref={(elem) => drawGraph(elem, open, 750, 750)} />
      </Modal>

    </>
  );
}

export default GeneMapperGraphs;
