import React, { useCallback } from 'react';
import { Stack, Container } from '@mui/material';
import Uploader from 'components/Uploader';
import StatusQueue from 'components/StatusQueue';
import styles from './dashboard.module.css';

function Dashboard() {
  return (
    <Stack
      direction="column"
      sx={{
        paddingTop: '100px',
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
