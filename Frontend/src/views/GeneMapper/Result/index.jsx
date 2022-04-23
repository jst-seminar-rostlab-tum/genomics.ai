import { ArrowBackIos, InfoOutlined } from '@mui/icons-material';
import ShareIcon from '@mui/icons-material/Share';
import {
  Button, Toolbar, Typography,
  Box,
  IconButton,
  Divider,
} from '@mui/material';
import React, { useCallback } from 'react';
/**
 * @param sidebarShown set true if sidebar is open, false otherwise
 */
function GeneMapperResultView({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '350px' : '100px'), [sidebarShown]);

  return (
    <Box sx={{ pl: paddingL, pr: '20px' }}>
      <Toolbar disableGutters>
        <Button startIcon={<ArrowBackIos fontSize="small" />} size="small" sx={{ mr: 2 }}>
          <Typography variant="caption" disableGutters>Back to projects</Typography>
        </Button>
        <Typography variant="h6">
          Projectname
        </Typography>
        <IconButton aria-label="learn more" size="small">
          <InfoOutlined fontSize="small" />
        </IconButton>
        <IconButton aria-label="share">
          <ShareIcon />
        </IconButton>
      </Toolbar>
      <Divider />
    </Box>
  );
}

export default GeneMapperResultView;
