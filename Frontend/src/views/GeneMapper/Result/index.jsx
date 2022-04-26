import { Box } from '@mui/material';
import GeneMapperCategories from 'components/GeneMapper/Categories';
import Sidepanel from 'components/GeneMapper/Sidepanel';
import GeneMapperResultHeader from 'components/GeneMapper/ResultHeader';
import React, { useCallback } from 'react';

/**
 * Shows the UMAP visualization for a given project.
 * @param sidebarShown set true if sidebar is open, false otherwise
 */
function GeneMapperResultView({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '350px' : '100px'), [sidebarShown]);

  console.log('render');

  return (
    <Box
      sx={{
        pl: paddingL,
        pr: '20px',
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <GeneMapperResultHeader projectName="Demo Projectname" />
      <Box
        sx={{
          display: 'flex', flexGrow: 1, justifyContent: 'space-between', alignItems: 'stretch',
        }}
      >
        <Sidepanel title="Categories">
          <GeneMapperCategories categories={['Cell type', 'Batch']} />
        </Sidepanel>
        <Box sx={{ flexGrow: 1 }} />
        <Sidepanel title="Graphs" collapseToRight />
      </Box>
    </Box>
  );
}

export default GeneMapperResultView;
