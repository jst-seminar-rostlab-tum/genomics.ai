import React, { useEffect, useState } from 'react';
import {
  List, Box, Collapse, IconButton,
} from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import StatusCard from './StatusCard/StatusCard';
import styles from './statusqueue.module.css';
import { filterDone, queryJobs } from './StatusQueueLogic';
import { JOB_QUEUE_UPDATE_INTERVAL } from '../../common/constants';

function StatusQueue() {
  const [expandList, setExpandList] = useState(true);
  const [jobs, setJobs] = useState([{
  _id: 1,
  status: 'DONE',
  }]);

  useEffect(() => {
    const updateJobs = () => queryJobs().then((newJobs) => setJobs(filterDone(newJobs)))
      .catch((ignored) => { console.log(ignored); });
    updateJobs();
    const intervalId = setInterval(updateJobs, JOB_QUEUE_UPDATE_INTERVAL);
    return (() => clearInterval(intervalId));
  }, [setJobs]);

  return (
    <div className={styles.StatusQueue}>
      <Box sx={{
        backgroundColor: 'rgba(0, 60, 255, 0.05)',
        width: '100%',
        padding: '10px',
        borderRadius: expandList ? '10px 10px 0px 0px' : '10px 10px 10px 10px',
        textAlign: 'center',
        fontWeight: 'bold',
        display: 'flex',
        justifyContent: 'center',
        fontSize: '8px',
      }}
      >
        <h1>Processing Queue</h1>
        <IconButton onClick={() => setExpandList(!expandList)}>
          {expandList ? <ExpandLessIcon /> : <ExpandMoreIcon />}
        </IconButton>
      </Box>

      <Box className={styles.flexContainer}>
        <Collapse in={expandList} sx={{ width: '100%', justifyContent: 'center', alignItems: 'center' }}>
          <List sx={{
            width: '100%', alignItems: 'center', justifyContent: 'center', flexDirection: 'column',
          }}
          >
            {jobs.map((job) => (
              <StatusCard key={job._id} id={job._id} status={job.status} log={job.fileName} sx={{ alignItems: 'center' }} />
            ))}
          </List>
        </Collapse>
      </Box>
    </div>
  );
}

export default StatusQueue;

/*
map jobs to list items containing the status cards
*/
