import React from 'react';
import { Box } from '@mui/material';
import Visualization from '../../Visualization/src/Visualization';
import styles from './visualizationpage.module.css';

function VisualizationPage({ id }) {

  //TODO: a new page is supposed to be created per completed project

  return (
    <Box sx={{
      display: 'flex',
      justifyContent: 'center',
      textAlign: 'center'
    }}
      className={styles.visualizationContainer}
    >
      <Visualization />
    </Box>
  );
}

export default VisualizationPage;
