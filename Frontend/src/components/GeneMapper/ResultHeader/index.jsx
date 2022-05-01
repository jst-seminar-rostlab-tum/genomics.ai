import { ArrowBackIos, InfoOutlined } from '@mui/icons-material';
import {
  Box, Divider, IconButton, Typography,
} from '@mui/material';
import CustomButton from 'components/CustomButton';
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
      <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
        <CustomButton type="tertiary" onClick={() => history.goBack()}>
          <ArrowBackIos sx={{ ml: -2 }} fontSize="small" />
          <Typography variant="caption">Back</Typography>
        </CustomButton>
        <Box sx={{ display: 'flex', alignItems: 'center' }}>
          <Typography variant="h6">
            {projectName}
          </Typography>
          <IconButton aria-label="learn more" size="small">
            <InfoOutlined fontSize="small" />
          </IconButton>
        </Box>
        {/* <Divider orientation="vertical" variant="middle" sx={{ m: 1 }} flexItem /> */}
        <ShareMenu projectName={projectName} url="https://www.genecruncher.com" />
      </Box>
      <Divider sx={{ mt: 1, mb: 1 }} />
    </>
  );
}

export default GeneMapperResultHeader;
