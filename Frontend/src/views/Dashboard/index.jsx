import { Container, Stack } from '@mui/material';
import React, { useCallback } from 'react';
import Uploader from 'components/Uploader';
import StatusQueue from 'components/StatusQueue';
import styles from './dashboard.module.css';

function Dashboard({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '100px' : '350px'), [sidebarShown]);
  return (
    <Stack
      direction="column"
      sx={{
        paddingTop: '10px',
        paddingLeft: paddingL,
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
