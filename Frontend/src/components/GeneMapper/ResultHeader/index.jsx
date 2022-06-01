import { ArrowBackIos, InfoOutlined } from '@mui/icons-material';
import DownloadIcon from '@mui/icons-material/Download';
import {
  Box, Divider, IconButton, Typography,
} from '@mui/material';
import CustomButton from 'components/CustomButton';
import { Modal } from 'components/Modal';
import React, { useEffect, useState } from 'react';
import { useHistory } from 'react-router-dom';
import AtlasService from 'shared/services/Atlas.service';
import ModelsService from 'shared/services/Models.service';
import ProjectInfo from '../ProjectInfo';
import ShareMenu from '../ShareMenu';

/**
 * Header shown on the result page, including navigation buttons, project info,
 * and sharing and download options
 * @param project Object containing project data
 */
function GeneMapperResultHeader({ project }) {
  const history = useHistory();
  const [showInfo, setShowInfo] = useState(false);
  const [atlas, setAtlas] = useState(null);
  const [model, setModel] = useState(null);

  useEffect(() => {
    AtlasService.getAtlasById(project.atlasId).then((data) => setAtlas(data));
    ModelsService.getModelById(project.modelId).then((data) => setModel(data));
  }, [project.atlasId, project.modelId]);

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
        <Box sx={{ p: 1 }}>
          <Typography variant="h4" sx={{ pb: 1 }}>{project.name}</Typography>
          <ProjectInfo project={project} atlas={atlas} model={model} />
        </Box>
      </Modal>
    </>
  );
}

export default GeneMapperResultHeader;
