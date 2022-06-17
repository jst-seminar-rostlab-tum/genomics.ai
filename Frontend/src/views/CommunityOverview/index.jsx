import React from 'react';
import HeaderView from 'components/general/HeaderView';
import styles from './communityOverview.module.css';
import TeamOverview from 'components/teams/overview/TeamOverview';
import InstitutionOverview from 'components/institutions/InstitutionOverview';
import { Divider } from '@mui/material';

function CommunityOverview({ sidebarShown }) {
  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="Community"
    >
      <div className={styles.content}>
        <TeamOverview />
        <br />
        <Divider />
        <br />
        <InstitutionOverview />
      </div>
    </HeaderView>
  );
}

export default CommunityOverview;
