import React, { useState } from 'react';
import {
  List, Box, Collapse, IconButton,
} from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import StatusCard from './StatusCard/StatusCard';
import styles from './statusqueue.module.css';

function StatusQueue() {
  const arr1 = [1, 2, 3, 4];
  const [expandList, setExpandList] = useState(true);

  return (
    <div className={styles.StatusQueue}>
      <div className={styles.listTitle}>
        <h1>Processing Queue</h1>
        <IconButton onClick={() => setExpandList(!expandList)}>
          {expandList ? <ExpandLessIcon /> : <ExpandMoreIcon />}
        </IconButton>
      </div>
      <Box className={styles.flexContainer}>
        <Collapse in={expandList} sx={{ width: '100%', justifyContent: 'center', alignItems: 'center' }}>
          <List sx={{ width: '100%', alignItems: 'center', justifyContent: 'center' }}>
            {arr1.map((index) => (
              <StatusCard id={index} sx={{ alignItems: 'center' }} />
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
