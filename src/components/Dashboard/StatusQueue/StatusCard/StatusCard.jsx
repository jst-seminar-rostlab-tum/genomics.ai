import React, { useState } from 'react';
import {
  Card,
  Box,
  Typography,
} from '@mui/material';
import {
  red,
  green,
  yellow,
  grey,
  blue,
} from '@mui/material/colors';
import CircleIcon from '@mui/icons-material/Circle';
import styles from './statuscard.module.css';

/*
The status card includes the status of the job.
Status types that the status card can include are:
pending - state during the upload of the file and during verification
processing - state during the processing of the file for the ML model
error - can occure in all steps including verification
completed - state after all stages are completed
unknown - the unknown state happens when the status can't be updated
*/

function StatusCard({ id }) {
  const [status, setStatus] = useState('unknown');
  const statusColor = {
    pending: yellow[600],
    processing: blue[300],
    error: red[300],
    completed: green[300],
    unknown: grey[500],
  };

  // setting status to known for testing
  const testStatusColor = statusColor.pending;
  //The status is not passed down as a prop, but received from backend
  setStatus('pending');

  return (
    <Box className={styles.cardContainer}>
      <Card variant="outlined">
        <div className={styles.headerContainer}>
          <Typography
            variant="h4"
            sx={{
              fontWeight: 'bold',
              fontSize: '16px',
              margin: '0',
              padding: '0',
              textAlign: 'center',
            }}
          >
            {`Job ${id}`}
          </Typography>
          <CircleIcon sx={{ color: testStatusColor, paddingLeft: '5px' }} />
        </div>
        <Typography variant="h4" sx={{ fontWeight: 'light', fontSize: '16px' }}>
          Status:
          <Typography
            variant="h4"
            sx={{
              fontSize: '16px',
              fontWeight: 'light',
              margin: '0',
              padding: '0',
              textAlign: 'center',
              display: 'inline',
              color: testStatusColor,
            }}
          >
            {` ${status}`}
          </Typography>
        </Typography>
      </Card>
    </Box>
  );
}

export default StatusCard;
