import { ArrowBackIos, InfoOutlined } from '@mui/icons-material';
import DownloadIcon from '@mui/icons-material/Download';
import {
  Box, Divider, IconButton, Typography,
} from '@mui/material';
import CustomButton from 'components/CustomButton';
import { Modal, ModalTitle } from 'components/Modal';
import React, { useState } from 'react';
import { useHistory } from 'react-router-dom';
import ShareMenu from '../ShareMenu';

/**
 *
 * @param projectName name of the visualized project
 */
function GeneMapperResultHeader({ project }) {
  const history = useHistory();
  const [showInfo, setShowInfo] = useState(false);

  return (
    <>
      <Box sx={{ display: 'flex', justifyContent: 'space-between' }}>
        <CustomButton type="tertiary" onClick={() => history.goBack()}>
          <ArrowBackIos sx={{ ml: -2 }} fontSize="small" />
          <Typography variant="caption">Back</Typography>
        </CustomButton>
        <Box sx={{ display: 'flex', alignItems: 'center' }}>
          <Typography variant="h6">
            {project.name}
          </Typography>
          <IconButton aria-label="learn more" size="small" onClick={() => setShowInfo(true)}>
            <InfoOutlined fontSize="small" />
          </IconButton>
        </Box>
        <Box sx={{ display: 'flex', alignItems: 'center' }}>
          <ShareMenu projectName={project.name} url={window.location} />
          <Divider orientation="vertical" flexItem variant="middle" />
          <IconButton href={project.location} download={`${project.name}.tsv`}>
            <DownloadIcon />
          </IconButton>
        </Box>
      </Box>
      <Divider sx={{ mt: 1, mb: 1 }} />
      <Modal isOpen={showInfo} setOpen={setShowInfo}>
        <ModalTitle>Test</ModalTitle>
        {Object.entries(project).map(([key, value]) => (<Typography key={key}>{`${key}: ${value}`}</Typography>))}
      </Modal>
    </>
  );
}

export default GeneMapperResultHeader;
