import { Box, Container } from '@mui/material';
import React from 'react';
import Uploader from './Uploader/Uploader';
import styles from './dashboard.module.css';

function Dashboard() {
  return (
    <Box
      className="flexContainer"
      sx={{
        display: 'flex',
        width: '100%',
        justifyContent: 'center',
      }}
    >
      <Container className={styles.sideBar}>
        <div><h1>Project Bar</h1></div>
      </Container>

      <Container className={styles.FileUpload}>
        <Uploader />
      </Container>

      <Container className={styles.fileQueue}>
        <h1>Queue of files</h1>
      </Container>
    </Box>
  );
}

//

export default Dashboard;
