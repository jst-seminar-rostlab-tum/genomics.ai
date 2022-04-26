import { ArrowBackIos, InfoOutlined } from '@mui/icons-material';
import {
  Button, Divider, IconButton, Toolbar, Typography,
} from '@mui/material';
import React from 'react';
import { useHistory } from 'react-router-dom';
import ShareMenu from '../ShareMenu';

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
        <Divider orientation="vertical" variant="middle" sx={{ m: 1 }} flexItem />
        <ShareMenu projectName={projectName} url="https://genecruncher.com" />
      </Toolbar>
      <Divider sx={{ mb: 1 }} />
    </>
  );
}

export default GeneMapperResultHeader;
