import { Container } from '@mui/material';
import React from 'react';
import HeaderView from 'components/HeaderView';
import ProfileImage from 'components/ProfileImage';
import styles from './institutionOverview.module.css';

function InstitutionOverview({ sidebarShown }) {
  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="My Institutions"
      replaceHeaderRight={<ProfileImage sizePixels={44} />}
      content={<Container className={styles.test}>Hello world</Container>}
    />
  );
}

export default InstitutionOverview;
