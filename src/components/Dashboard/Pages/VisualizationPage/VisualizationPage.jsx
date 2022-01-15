import React from 'react';
import { Box } from '@mui/material';
import Visualization from '../../Visualization/src/Visualization';
import styles from './visualizationpage.module.css';

function VisualizationPage({ id }) {
  /*
  TODO: a new page is supposed to be created per completed project
  the presented page is dependent on the result
  See how to implement the id inside of the visualization project result
  */
  return (
    <Box
      sx={{
        display: 'flex',
        justifyContent: 'center',
        textAlign: 'center',
      }}
      className={styles.visualizationContainer}
    >
      <Visualization id={id} />
    </Box>
  );
}

export default VisualizationPage;
