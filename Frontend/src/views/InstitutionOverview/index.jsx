import { Container } from '@mui/material';
import React from 'react';
import HeaderView from 'components/HeaderView';
import styles from './institutionOverview.module.css';

function InstitutionOverview({ sidebarShown }) {
  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="My Institutions"
    >
      <Container className={styles.test}>Hello world</Container>
    </HeaderView>
  );
}

export default InstitutionOverview;
