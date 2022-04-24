import { Box } from '@mui/material';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import React, { useCallback } from 'react';

/**
 * Shows the UMAP visualization for a given project.
 * @param sidebarShown set true if sidebar is open, false otherwise
 */
function GeneMapperResultView({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '350px' : '100px'), [sidebarShown]);

  return (
    <Box sx={{ pl: paddingL, pr: '20px' }}>
      <GeneMapperResultHeader projectName="Demo Projectname" />
    </Box>
  );
}

export default GeneMapperResultView;
