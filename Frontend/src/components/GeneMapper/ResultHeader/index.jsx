import { ArrowBackIos, InfoOutlined } from '@mui/icons-material';
import ShareIcon from '@mui/icons-material/Share';
import {
  Button, Divider, IconButton, Toolbar, Typography,
} from '@mui/material';
import React from 'react';

import { useHistory } from 'react-router-dom';

/**
 *
 * @param projectName name of the visualized project
 */
function GeneMapperResultHeader({ projectName }) {
  const history = useHistory();

  return (
    <>
      <Toolbar disableGutters>
        <Button startIcon={<ArrowBackIos fontSize="small" />} size="small" sx={{ mr: 2 }} onClick={() => history.goBack()}>
          <Typography variant="caption" disableGutters>Back to projects</Typography>
        </Button>
        <Typography variant="h6">
          {projectName}
        </Typography>
        <IconButton aria-label="learn more" size="small">
          <InfoOutlined fontSize="small" />
        </IconButton>
        <IconButton aria-label="share">
          <ShareIcon />
        </IconButton>
      </Toolbar>
      <Divider sx={{ mb: 1 }} />
    </>
  );
}

export default GeneMapperResultHeader;
