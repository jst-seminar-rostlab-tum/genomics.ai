import React from 'react';
import { Typography } from '@mui/material';

function ProjectInfo({ project, atlas, model }) {
  return (
    <>
      <Typography>{`Atlas: ${atlas?.name}`}</Typography>
      <Typography>{`Model: ${model?.name}`}</Typography>
      <Typography>{`Dataset: ${project?.fileName}`}</Typography>
    </>
  );
}

export default ProjectInfo;
