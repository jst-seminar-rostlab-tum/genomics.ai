import React, { useState } from 'react';
import {
  List, Box, Collapse, IconButton,
} from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import StatusCard from './StatusCard/StatusCard';
import styles from './statusqueue.module.css';

function StatusQueue() {
  const arr1 = [1, 2, 3, 4, 5, 6];
  const [expandList, setExpandList] = useState(true);

  return (
    <div className={styles.StatusQueue}>
      <Box sx={{
        backgroundColor: 'rgba(0, 60, 255, 0.05)',
        width: '100%',
        padding: '10px',
        borderRadius: expandList ? '10px 10px 0px 0px' : '10px 10px 10px 10px',
        textAlign: 'center',
        fontEeight: 'bold',
        display: 'flex',
        justifyContent: 'center',
        fontSize: '8px'
      }}>
        <h1>Processing Queue</h1>
        <IconButton onClick={() => setExpandList(!expandList)}>
          {expandList ? <ExpandLessIcon /> : <ExpandMoreIcon />}
        </IconButton>
      </Box>



      <Box className={styles.flexContainer}>
        <Collapse in={expandList} sx={{ width: '100%', justifyContent: 'center', alignItems: 'center' }}>
          <List sx={{
            width: '100%', alignItems: 'center', justifyContent: 'center', flexDirection: 'column'
          }}
          >
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
