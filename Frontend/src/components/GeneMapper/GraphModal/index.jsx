import { Box } from '@mui/material';
import React from 'react';

function GraphModal({ draw }) {
  return (
    <Box sx={{ backgroundColor: 'white' }}>
      <Box ref={(elem) => draw(elem, 500, 500)} />
    </Box>
  );
}

export default GraphModal;
