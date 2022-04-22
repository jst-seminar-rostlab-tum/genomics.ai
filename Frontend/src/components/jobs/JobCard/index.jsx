import React from 'react';
import Circle from '@mui/icons-material/Circle';
import Button from '@mui/material/Button';
import ListCard from 'components/general/ListCard';
import { jobStatusColors } from 'shared/utils/common/constants';
// import styles from './teamCard.module.css';

function JobCard({ job }) {
  const navigateToJob = () => window.open(job.location, '_blank');

  const {
    id, status,
  } = job;
  return (
    <div onClick={navigateToJob} onKeyPress={navigateToJob} role="button" tabIndex={-1}>
      <ListCard
        imageComponent={(
          <Circle sx={{ color: jobStatusColors[status], width: '32px', height: '32px' }} />
        )}
        title={`Job ${id.substring(id.length - 4)}`}
        trailing={(
          <Button>See Results</Button>
        )}
      />
    </div>
  );
}

export default JobCard;
