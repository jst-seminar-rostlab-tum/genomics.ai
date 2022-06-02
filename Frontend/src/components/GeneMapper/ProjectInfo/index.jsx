import React, { useState } from 'react';
import { IconButton, Typography } from '@mui/material';
import AtlasInfo from '../AtlasInfo';
import ModelInfo from '../ModelInfo';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
/**
 * Displays information about the given project
 * @param project Object containing project data
 * @param atlas Object containing atlas data
 * @param model Object containing model data
 */
function ProjectInfo({ project, atlas, model }) {
  const [atlasInfoOpen, setAtlasInfoOpen] = useState(false);
  const [modelInfoOpen, setModelInfoOpen] = useState(false);

  return (
    <>
      <Typography>
        {`Atlas: ${atlas?.name}`}
        <IconButton size="small" onClick={() => setAtlasInfoOpen(true)}>
          <InfoOutlinedIcon fontSize="small" />
        </IconButton>
      </Typography>
      <Typography>
        {`Model: ${model?.name}`}
        <IconButton size="small" onClick={() => setModelInfoOpen(true)}>
          <InfoOutlinedIcon fontSize="small" />
        </IconButton>
      </Typography>
      <Typography>{`Dataset: ${project?.fileName}`}</Typography>
      {atlas
      && <AtlasInfo id={atlas._id} open={atlasInfoOpen} setOpen={setAtlasInfoOpen} />}
      {model
      && <ModelInfo id={model._id} open={modelInfoOpen} setOpen={setModelInfoOpen} />}
    </>
  );
}

export default ProjectInfo;
