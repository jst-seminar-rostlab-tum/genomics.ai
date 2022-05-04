import React from 'react';
import Circle from '@mui/icons-material/Circle';
import ListCard from 'components/general/ListCard';
import { jobStatusColors } from 'shared/utils/common/constants';
import CustomButton from 'components/CustomButton';
// import styles from './teamCard.module.css';

function JobCard({ job }) {
  const navigateToJob = () => window.open(job.location, '_blank');

  const {
    _id, status,
  } = job;
  return (
    <div onClick={navigateToJob} onKeyPress={navigateToJob} role="button" tabIndex={-1}>
      <ListCard
        imageComponent={(
          <Circle sx={{ color: jobStatusColors[status], width: '32px', height: '32px' }} />
        )}
        title={`Job ${_id}`}
        trailing={(
          <CustomButton type="primary">See Results</CustomButton>
        )}
      />
    </div>
  );
}

export default JobCard;
