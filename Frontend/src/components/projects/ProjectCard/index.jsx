import React from 'react';
import Circle from '@mui/icons-material/Circle';
import ListCard from 'components/general/ListCard';
import { jobStatusColors } from 'shared/utils/common/constants';
import CustomButton from 'components/CustomButton';
// import styles from './teamCard.module.css';

function ProjectCard({ project }) {
  const navigateToProject = () => window.open(project.location, '_blank');

  const {
    _id, status,
  } = project;
  return (
    <div onClick={navigateToProject} onKeyPress={navigateToProject} role="button" tabIndex={-1}>
      <ListCard
        imageComponent={(
          <Circle sx={{ color: jobStatusColors[status], width: '32px', height: '32px' }} />
        )}
        title={`Project ${_id}`}
        trailing={(
          <CustomButton type="primary">See Results</CustomButton>
        )}
      />
    </div>
  );
}

export default ProjectCard;
