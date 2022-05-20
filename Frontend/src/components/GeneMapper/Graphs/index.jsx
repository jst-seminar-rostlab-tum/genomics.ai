import React, { useEffect, useRef } from 'react';
import { Box } from '@mui/material';

function GeneMapperGraphs({ graphs, drawGraph }) {
  const graphContainers = useRef({});

  useEffect(() => {
    Object.entries(graphContainers.current).forEach(([id, elem]) => {
      drawGraph(elem, id, 270, 270);
    });
  }, [graphContainers.current]);

  return (
    <Box sx={{ height: '100%' }}>
      {graphs.map((id) => (
        <Box key={id} ref={(elem) => { graphContainers.current[id] = elem; }} />
      ))}

    </Box>
  );
}

export default GeneMapperGraphs;
