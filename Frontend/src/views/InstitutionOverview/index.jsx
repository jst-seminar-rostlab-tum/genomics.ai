import { Container, Stack } from '@mui/material';
import React, { useCallback } from 'react';
import Header from 'components/Header';
import styles from './institutionOverview.module.css';

function InstitutionOverview({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '100px' : '350px'), [sidebarShown]);
  return (
    <Stack
      direction="column"
      sx={{
        paddingLeft: paddingL,
        display: 'flex',
        width: '100%',
        justifyContent: 'center',
      }}
    >
      <Header title="My Institutions" />

      <Stack
        className="flexContainer"
        direction="row"
      >
        <Container className={styles.test}>
          Hello
        </Container>
      </Stack>
    </Stack>
  );
}

export default InstitutionOverview;
