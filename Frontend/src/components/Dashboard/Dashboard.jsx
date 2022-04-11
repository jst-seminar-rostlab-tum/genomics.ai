import { Container, Stack } from '@mui/material';
import React from 'react';
import Uploader from './Uploader/Uploader';
import StatusQueue from './StatusQueue/StatusQueue';
import styles from './dashboard.module.css';

function Dashboard() {
  return (
    <Stack
      direction="column"
      sx={{
        paddingTop: '40px',
        paddingLeft: '100px',
        display: 'flex',
        width: '100%',
        justifyContent: 'center',
      }}
    >
      <div className={styles.title}>
        <h1>Dashboard</h1>
      </div>

      <Stack
        className="flexContainer"
        direction="row"
      >
        <Container className={styles.fileUpload}>
          <Uploader />
        </Container>

        <Container className={styles.fileQueue}>
          <StatusQueue />
        </Container>
      </Stack>
    </Stack>
  );
}

export default Dashboard;
