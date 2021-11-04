import React, {useState} from 'react';
import { List, ListItem, Box, Collapse, IconButton } from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import StatusCard from './StatusCard/StatusCard';
import styles from './statusqueue.module.css';

function StatusQueue() {
  let arr1 = [1, 2, 3];
  const [expandList, setExpandList] = useState(false);

  return (
    <div>
      <div className={styles.listTitle}>
        <h1>Processing Queue</h1>
        <IconButton onClick={() => setExpandList(!expandList)}>
          <ExpandMoreIcon />
        </IconButton>
      </div>
      <Box className={styles.flexContainer}>
        <Collapse in={expandList}>
          <List>
            {arr1.map((index) => (
              <StatusCard 
              index={index}
              ></StatusCard>
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